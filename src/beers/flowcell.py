import gzip
import itertools
import numpy as np
import re
import pysam
import math
import os
from dataclasses import dataclass, field
from typing import Generator

from beers.cluster_packet import ClusterPacket
from beers.beers_exception import BeersException
from beers.cluster import Cluster


@dataclass
class FlowcellLane:
    """
    These objects partly compose the flowcell object, representing the lanes that the current BEERS run is using.  The
    object keeps track of all used coordinates to prevent duplication.

    Attributes
    ---------
    lane
        Lane number
    consumed_coordinates
        set of coordinates already used
    """

    lane: int
    consumed_coordinates: set = field(default_factory = set)

coords_match_pattern = re.compile(r'^.*:(\w+):(\d+):(\d+):(\d+):(\d+)$')
@dataclass
class Flowcell:
    """
    Models a flowcell. MoleculePackets can be loaded and are converted to ClusterPackets with assigned coordinates
    for all the clusters.

    A flowcell's geometry (i.e., coordinate ranges for lane, tile, x, y) are either provided via the configuration
    The flowcell retention rate is specified as a percentage in the configuration file.  The user may restrict which
    lanes of the flowcell to use via the configuration. Otherwise, all the lanes described by the flowcell's geometry
    are used.  The individual flowcell lanes to be used are instantiated as FlowcellLane objects, each of which keeps
    track of the coordinates used for its given lane.  Each lane has its own coordinate generator to produce new
    coordinates for each newly retained molecule.
    """

    parameters: dict
    rng: np.random.Generator
    min_coords: dict = field(default_factory = lambda: {"lane": 0, "tile": 0, "x": 0, "y": 0})
    max_coords: dict = field(default_factory = lambda: {"lane": 10_000, "tile": 10_000, "x": 10_000, "y": 10_000})
    available_lanes: list[int] = field(default_factory = list)
    lanes_to_use: list[int] = field(default_factory = list)
    flowcell_lanes: dict[str, FlowcellLane] = field(default_factory = dict)
    coordinate_generators: dict[str, Generator] = field(default_factory = dict)

    def __post_init__(self):
        """
        """
        self.set_flowcell_coordinate_ranges()
        self.available_lanes = list(range(self.min_coords['lane'], self.max_coords['lane'] + 1))
        self.lanes_to_use = self.parameters["lanes_to_use"] or self.available_lanes

        if self.parameters['coordinate_strategy'] == 'random':
            # Generate coordinates with replacement, i.e. may not be distinct
            # this is the fastest/simplest method
            generate_coordinates = self.generate_coordinates_random
        elif self.parameters['coordinate_strategy'] == 'random_distinct':
            # Generate coordinates without replacement, i.e. are all distinct
            # Uses much more memory if sequencing millions of reads
            generate_coordinates = self.generate_coordinates_distinct

        for lane in self.lanes_to_use:
            self.flowcell_lanes[lane] = FlowcellLane(lane)
            self.coordinate_generators[lane] = generate_coordinates(lane)

    def validate(self):
        """
        Insures the the configuration file parameters relating to the flowcell are appropriate.  Otherwise one or more
        error messages are crafted.
        :return: A tuple indicating whether or not the flowcell configuration is valid and a diagnostic message for
        the user.
        """
        valid = True
        msg = ""
        if not set(self.lanes_to_use).issubset(set(self.available_lanes)):
            valid = False
            msg += f"The flowcell lanes to use {self.lanes_to_use} must be a subset of the available lanes" \
                   f" {self.available_lanes}.\n"
        if self.parameters['coordinate_strategy'] not in ['random', 'random_distinct']:
            valid = False
            msg += f"Flowcell coordinate_strategy must be one of 'random' or 'random_distinct'"
        return valid, msg

    def convert_molecule_pkt_to_cluster_pkt(self, molecule_packet):
        """
        Takes the molecule packet containing only retained molecules and evenly apportions those molecules across the
        flowcell lanes available for use (i.e., first group is attached to the first lane to use, the seocnd group to
        the second lane to use and so on).  Each molecule is then used to create a cluster, identified by the lane
        and the tile,x, and y coordinates as provided by the lane's coordinate generator.  The clusters are rolled
        up into a cluster packet, which again picks up the molecule packet's sample information.  So the cluster
        packet, like the molecule packet from which it was derived, also applies to one sample only.
        :param molecule_packet: The molecule packet now containing only those molecule retained by the flowcell.
        :return: A cluster packet containing the retained molecules as cluster objects with a lane and set of
        coordinates identified.
        """
        cluster_packet_id = molecule_packet.molecule_packet_id
        clusters = []
        molecules_per_lane = len(molecule_packet.molecules)//len(self.lanes_to_use)
        lane_index = 0
        lane = self.lanes_to_use[lane_index]
        for counter, molecule in enumerate(molecule_packet.molecules):
            cluster_id = Cluster.next_cluster_id
            if (counter + 1) % molecules_per_lane == 0 and lane_index + 1 < len(self.lanes_to_use):
                lane_index += 1
                lane = self.lanes_to_use[lane_index]
            clusters.append(Cluster(cluster_id, molecule, lane, next(self.coordinate_generators[lane])))
            Cluster.next_cluster_id += 1
        print(f"Assigned flowcell coordinates to {counter + 1} clusters.")
        return ClusterPacket(cluster_packet_id, molecule_packet.sample, clusters)

    def load_flowcell(self, molecule_packet):
        """
        The entry point into the flowcell object.  Aside from the instantiation and validation, all other bound methods
        are helper methods.  The portion of the molecules of the molecule packet retained by the flowcell replaces the
        original molecule collection and that diminished molecule packet is converted into a cluster packet, where
        each molecule is converted into a cluster having a unique set of flowcell coordinates.
        :param molecule_packet: incoming molecule packet (from library prep or a simulated library prep output)
        :return: cluster packet contained those retained molecules affixed to the flowcell via unique coordinates.
        """
        cluster_packet = self.convert_molecule_pkt_to_cluster_pkt(molecule_packet)
        return cluster_packet

    def generate_coordinates_distinct(self, lane):
        """
        A Python generator that randomly selects tile, x, and y coordinates and uses them to create a  coordinates
        tuple (tile, x, y).  The new lane coordinates object is compared with all the consumed coordinates in the given lane.  If
        the object is in the consumed coordinates, it is discarded and the random selection repeated.  Once a selection
        is found to be unique, it is added to the consumed coordinates for the given lane and the object is returned.
        If no set of unique coordinates can be found in 100 attempts, an exception is thrown.
        :param lane: the lane for which a unique combination of tile, x, and y cooredinates is to be selected
        :return: a unique combination of coordinates for the lane given.
        """
        ctr = 0
        flowcell_lane = self.flowcell_lanes[lane]
        consumed_coordinates = flowcell_lane.consumed_coordinates
        while True:
            ctr += 1
            x = self.rng.integers(self.min_coords['x'], self.max_coords['x'] + 1)
            y = self.rng.integers(self.min_coords['y'], self.max_coords['y'] + 1)
            tile = self.rng.integers(self.min_coords['tile'], self.max_coords['tile'] + 1)
            coord = (tile, x, y)
            if coord not in consumed_coordinates:
                ctr = 0
                consumed_coordinates.add(coord)
                yield coord
            if ctr >= 100:
                raise BeersException("Unable to find unused flowcell coordinates after 100 attempts.")

    def  generate_coordinates_random(self, lane):
        """
        A Python generator that randomly selects tile, x, and y coordinates and uses them to create a  coordinates
        tuple (tile, x, y). May generate the same coordinate location multiple times, so not appropriate for simulations
        that are intended for uses that depend upon coordinates. Lower memory usage than generate_coordinates_distinct.
        :param lane: the lane for which a unique combination of tile, x, and y cooredinates is to be selected
        :return: a unique combination of coordinates for the lane given.
        """
        flowcell_lane = self.flowcell_lanes[lane]
        while True:
            x = self.rng.integers(self.min_coords['x'], self.max_coords['x'] + 1)
            y = self.rng.integers(self.min_coords['y'], self.max_coords['y'] + 1)
            tile = self.rng.integers(self.min_coords['tile'], self.max_coords['tile'] + 1)
            coord = (tile, x, y)
            yield coord

    def set_flowcell_coordinate_ranges(self):
        """
        Determines the minimum and maximum range for each flowcell coordinate, including the lane.  If the ranges are
        supplied via the configuration file, those ranges (min and max) are rendered as dictionaries.  Otherwise, each
        input FASTQ file is applied to obtain the same ranges.
        """
        geometry = self.parameters['flowcell_geometry']
        self.min_coords = {"lane": geometry['min_lane'], "tile": geometry['min_tile'],
                           "x": geometry['min_x'], "y": geometry['min_y']}
        self.max_coords = {"lane": geometry['max_lane'], "tile": geometry['max_tile'],
                           "x": geometry['max_x'], "y": geometry['max_y']}
