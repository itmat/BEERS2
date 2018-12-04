from collections import namedtuple


class FlowcellLane:

    def __init__(self, lane):
        self.consumed_coordinates = []
        self.lane = lane


class LaneCoordinates:

    def __init__(self, lane, tile, x, y):
        self.lane = lane
        self.tile = tile
        self.x = x
        self.y = y

    def __eq__(self, other):
        return self.lane == other.lane and self.tile == other.tile and self.x == other.x and self.y == other.y

    def __ne__(self, other):
        return not self.__eq__(other)

    def __gt__(self, other):
        if self.lane != other.lane:
            return self.lane > other.lane
        if self.tile != other.tile:
            return self.tile > other.tile
        if self.x != other.x:
            return self.x > other.x
        return self.y > other.y

    def __ge__(self, other):
        return self.__gt__(other) or self.__eq__(other)

    def __lt__(self, other):
        return not self.__gt__(other) and self.__ne__(other)

    def __le__(self, other):
        return not self.__gt__(other)

    def __str__(self):
        return f"lane: {self.lane}, tile: {self.tile}, x: {self.x}, y: {self.y}"

    def serialize(self):
        return f"{self.lane}\t{self.tile}\t{self.x}\t{self.y}"

    @staticmethod
    def deserialize(data):
        data = data[1:] if (data.startswith("#")) else data
        lane, tile, x, y = data.rstrip().split("\t")
        return LaneCoordinates(lane, tile, x, y)
