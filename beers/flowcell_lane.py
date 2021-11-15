from collections import namedtuple
import functools


class FlowcellLane:
    """
    These objects partly compose the flowcell object, representing the lanes that the current BEERS run is using.  The
    object keeps track of all used coordinates to prevent duplication.
    """

    def __init__(self, lane):
        """
        Instantiates a flowcell lane object for the given lane starting with no used coordinates.
        :param lane: the lane to which this object applies.
        """
        self.consumed_coordinates = set()
        self.lane = lane
