#!/usr/bin/env python3
import sys
from collections import Counter
from ntlink_patch_gaps import parse_minimizers
class Mapping:
    def __init__(self, read, ctg, num_anchors, mapping_str):
        self.read = read
        self.ctg = ctg
        self.num_anchors = num_anchors
        self.mappings = parse_minimizers(mapping_str)

    def __eq__(self, other):
        return self.read == other.read and self.ctg == other.ctg and \
               self.num_anchors == other.num_anchors and \
               self.mappings == other.mappings

    def equal_with_fewer_mappings(self, other):
        if self.read != other.read or self.ctg != other.ctg:
            return False
        if not set(self.mappings).issubset(set(other.mappings)):
            return False
        if self.num_anchors > other.num_anchors:
            return False
        return True

    def equal_with_more_mappings(self, other):
        if self.read != other.read or self.ctg != other.ctg:
            return False
        if not set(other.mappings).issubset(set(self.mappings)):
            return False
        if self.num_anchors < other.num_anchors:
            return False
        return True

def read_mappings(filename):
    mappings = {}
    with open(filename, "r") as f:
        for line in f:
            read, ctg, num_anchors, mappings_str = line.strip().split("\t")
            if read not in mappings:
                mappings[read] = []
            mappings[read].append(Mapping(read, ctg, num_anchors, mappings_str))
    return mappings

def assess_liftover_mappings(mappings, liftover_mapping_filename):
    c = Counter()
    with open(liftover_mapping_filename, "r") as f:
        for line in f:
            read, ctg, num_anchors, mappings_str = line.strip().split("\t")
            liftover_mapping = Mapping(read, ctg, num_anchors, mappings_str)
            if read not in mappings:
                c["not_in_full_rounds"] += 1
                continue
            has_full_equal = False
            has_addtl_mappings = False
            has_fewer_mappings = False
            for mapping in mappings[read]:
                if mapping == liftover_mapping:
                    has_full_equal = True
                if liftover_mapping.equal_with_more_mappings(mapping):
                    has_addtl_mappings = True
                if liftover_mapping.equal_with_fewer_mappings(mapping):
                    has_fewer_mappings = True
            if has_full_equal:
                c["full_equal"] += 1
            elif has_addtl_mappings:
                c["addtl_mappings"] += 1
            elif has_fewer_mappings:
                c["fewer_mappings"] += 1
            else:
                c["other"] += 1

    for category, count in c.items():
        print(f"{category}\t{count}")

def main():
    if len(sys.argv[1:]) != 2:
        print(f"Usage: {sys.argv[0]} <full rounds mapping> <liftover mappings>")
        sys.exit(1)

    full_rounds_mapping_filename = sys.argv[1]
    liftover_mapping_filename = sys.argv[2]

    mappings = read_mappings(full_rounds_mapping_filename)
    assess_liftover_mappings(mappings, liftover_mapping_filename)


if __name__ == "__main__":
    main()
