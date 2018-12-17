from collections import namedtuple


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
        self.consumed_coordinates = []
        self.lane = lane


class LaneCoordinates:
    """
    This object holds a set of tile, x, and y coordinates and is primarily used to allow for coordinate comparisons
    to insure again duplicate coordinates.  Since the object becomes part of a cluster and since a cluster via
    cluster packet is de/serialized, this object must also be serializable.
    """

    def __init__(self, tile, x, y):
        self.tile = tile
        self.x = x
        self.y = y

    def __eq__(self, other):
        """
        Objects are equivalent if the tile, x, and y coordinates are each equivalent.
        :param other: The LaneCoordinates object against which to compare.
        :return: True if equivalent and false otherwise
        """
        return self.tile == other.tile and self.x == other.x and self.y == other.y

    def __ne__(self, other):
        return not self.__eq__(other)

    def __gt__(self, other):
        """
        Familiar method for establishing which object is greater, assuming tiles are most significant than
        x's are more significant than y's.  Not needed to establish whether a pair of lane coordinate object is
        identical.  But is used to order resulting FASTQ files by coordinates. All other comparison methods are
        derived from this method and the eqnivalence method.
        :param other: The LaneCoordinates object against which to compare.
        :return: True if this object is greater and false otherwise.
        """
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
        return f"tile: {self.tile}, x: {self.x}, y: {self.y}"

    def serialize(self):
        """
        Simple serialization string
        :return: tab delimited string of coordinate set elements
        """
        return f"{self.tile}\t{self.x}\t{self.y}"

    @staticmethod
    def deserialize(data):
        """
        Removes any leading hash from the data string, if present, unpacks the remaining string into tile, x, and
        y and creates a new LaneCoordinates object from those elements.
        :param data: The string representation of the lane coordinates object
        :return: The populated LaneCoordinates object
        """
        data = data[1:] if (data.startswith("#")) else data
        tile, x, y = data.rstrip().split("\t")
        return LaneCoordinates(tile, x, y)
