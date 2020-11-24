# Defining helper classes
class PathNode:
    "Defines a node in a path of contig regions"
    def __init__(self, contig, ori, gap_size=None):
        self.contig = contig
        self.ori = ori
        self.gap_size = gap_size

    def set_gap_size(self, gap_size):
        "Set the gap size of the path node"
        self.gap_size = gap_size

    def get_ori_contig(self):
        return "{}{}".format(self.contig, self.ori)

    def get_gap(self):
        if self.gap_size is not None:
            return "{}N".format(self.gap_size)
        return None

    def __str__(self):
        return "{}{}".format(self.contig, self.ori)
