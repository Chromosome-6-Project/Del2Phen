#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Chromosome 6 Project - Utilities.

@author: T.D. Medina
"""

REFERENCE_CHR = [str(i) for i in range(1, 23)] + ["X", "Y"]


def jaccard(set1, set2):
    """Calculate Jaccard Index of two sets."""
    if isinstance(set1, (int, float)) and isinstance(set2, (int, float)):
        if set1 == 0 or set2 == 0:
            return 0
        return min([set1, set2]) / max([set1, set2])
    set1 = set(set1)
    set2 = set(set2)
    union = set1 | set2
    if len(union) == 0:
        return 0, set()
    intersect = set1 & set2
    jaccard_index = len(intersect) / len(union)
    return jaccard_index, intersect


def cnv_overlap(cnv1, cnv2):
    """Check if two CNVs overlap each other."""
    if cnv1.chromosome != cnv2.chromosome:
        return False
    return overlap(cnv1.range, cnv2.range)


def overlap(range1: range, range2: range):
    """Test if two ranges intersect."""
    # Needs to be <, not <=, because of range half-open notation.
    # e.g. range(1,3) should not overlap range(3, 5) because
    # range(1,3) is [1,2] and range(3,5) is [3,4].
    if range1.start < range2.stop and range2.start < range1.stop:
        return True
    return False


def get_range_intersect(range1: range, range2: range):
    """Make range from the intersect of two ranges."""
    if not overlap(range1, range2):
        return range(0)
    start = max(range1.start, range2.start)
    stop = min(range1.stop, range2.stop)
    return range(start, stop)


def merge_ranges(range1: range, range2: range):
    """Make range from the union of two ranges."""
    if not overlap(range1, range2):
        raise ValueError("Ranges do not overlap.")
    merged = range(min(range1.start, range1.start),
                   max(range1.stop, range2.stop))
    return merged


def merge_range_list(ranges):
    """Merge all overlapping ranges in a list of ranges."""
    ranges_copy = sorted(ranges.copy(), key=lambda x: x.stop)
    ranges_copy = sorted(ranges_copy, key=lambda x: x.start)
    merged = []

    while ranges_copy:
        range1 = ranges_copy[0]
        del ranges_copy[0]
        merges = []
        for i, range2 in enumerate(ranges_copy):
            if overlap(range1, range2):
                range1 = merge_ranges(range1, range2)
                merges.append(i)
        merged.append(range1)
        for i in reversed(merges):
            del ranges_copy[i]
    merged = sorted(merged, key=lambda x: x.start)

    return merged


def length_of_range_intersects(ranges):
    """Calculate total length of all range intersects in list of ranges."""
    ranges_copy = sorted(ranges.copy(), key=lambda x: x.stop)
    ranges_copy = sorted(ranges_copy, key=lambda x: x.start)
    intersects = []

    while ranges_copy:
        range1 = ranges_copy[0]
        del ranges_copy[0]
        for range2 in ranges_copy:
            if not overlap(range1, range2):
                continue
            intersects.append(get_range_intersect(range1, range2))
    intersects = merge_range_list(intersects)
    total = sum([len(x) for x in intersects])
    return total
