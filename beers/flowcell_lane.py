from collections import namedtuple

LaneCoordinates = namedtuple('Coordinates', ['tile', 'x', 'y'])

class FlowcellLane:

    def __init__(self, lane):
        self.consumed_coordinates = []
        self.lane = lane