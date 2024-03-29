'''
Helper functions for PAF output for ntLink mappings
'''

__author__ = "Lauren Coombe @lcoombe"

import io

def is_consistent(ctg_positions: list, increasing: bool, i1: int, i2: int, duplicate_positions: set) -> bool:
    "Given the expected monotonicity and a list with contig positions of mappings," \
    "determine whether the given indices are consistent"
    if ctg_positions[i1].ctg_pos in duplicate_positions or ctg_positions[i2].ctg_pos in duplicate_positions:
        return True  # Err on the side of caution when checking consistency including dups, just filter
    if increasing:
        return ctg_positions[i1].read_pos <= ctg_positions[i2].read_pos
    return ctg_positions[i1].read_pos >= ctg_positions[i2].read_pos

def break_mapping_blocks(sorted_ctg_pos: list, breaks: set, filters: set) -> list:
    "Break the mapping blocks at detected locations"
    return_mapping_blocks = []
    current_mapping_block = []
    for i, mapping in enumerate(sorted_ctg_pos):
        if i in filters:
            continue
        if i in breaks:
            return_mapping_blocks.append(current_mapping_block)
            current_mapping_block = [mapping]
        else:
            current_mapping_block.append(mapping)
    return_mapping_blocks.append(current_mapping_block)

    return return_mapping_blocks

def filter_and_break_mapping_blocks(transitions: list, sorted_ctg_pos: list, duplicate_positions: set,
                                    increasing=True) -> list:
    "Go through transitions in mapping block, find places to break block (if needed), and return broken blocks"
    breaks = set()
    filters = set()
    for i, transition in enumerate(transitions):
        if not transition:
            if sorted_ctg_pos[i].ctg_pos in duplicate_positions or sorted_ctg_pos[i + 1].ctg_pos in duplicate_positions:
                continue  # Just skip when transitions include duplicate contig positions
            if i + 2 >= len(transitions):
                # This is an end minimizer that's an issue. Remove it
                breaks.add(i + 1)
            elif is_consistent(sorted_ctg_pos, increasing, i, i + 2, duplicate_positions):
                # This is a single issue minimizer. Remove it
                filters.add(i + 1)
            elif i > 0 and is_consistent(sorted_ctg_pos, increasing, i - 1, i + 1, duplicate_positions):
                # This is a single issue minimizer. Remove it
                filters.add(i)
            else:
                # This is a larger segment problem or problem minimizer at the beginning. Break the mapping block
                breaks.add(i + 1)

    if not breaks and not filters:
        return [sorted_ctg_pos]
    return break_mapping_blocks(sorted_ctg_pos, breaks, filters)

def get_mapped_blocks(sorted_ctg_pos: list, min_consistent=0.75) -> list:
    "Go through a list of transitions, deciding whether to filter or cut the mapping block accordingly"
    ctg_positions = set()
    dup_positions = set()
    transitions_incr, transitions_decr = [], []
    all_transitions_incr, all_transitions_decr = True, True

    for i, j in zip(sorted_ctg_pos, sorted_ctg_pos[1:]):
        incr = i.read_pos <= j.read_pos
        transitions_incr.append(incr)
        all_transitions_incr = False if not incr else all_transitions_incr

        decr = i.read_pos >= j.read_pos
        transitions_decr.append(decr)
        all_transitions_decr = False if not decr else all_transitions_decr

        if i.ctg_pos in ctg_positions:
            dup_positions.add(i.ctg_pos)
        else:
            ctg_positions.add(i.ctg_pos)
    if sorted_ctg_pos[-1].ctg_pos in ctg_positions:
        dup_positions.add(sorted_ctg_pos[-1].ctg_pos)

    if all_transitions_incr or all_transitions_decr:
        return [sorted_ctg_pos]

    transition_counts_incr = transitions_incr.count(True)
    if (transition_counts_incr / len(transitions_incr)) >= min_consistent:
        return filter_and_break_mapping_blocks(transitions_incr, sorted_ctg_pos, dup_positions,
                                               increasing=True)
    if ((len(transitions_incr) - transition_counts_incr) / len(transitions_incr)) >= min_consistent:
        return filter_and_break_mapping_blocks(transitions_decr, sorted_ctg_pos, dup_positions,
                                               increasing=False)
    return []

def check_if_must_check_mapped_blocks(unsorted_hits: list, sorted_hits: list) -> bool:
    "Use sorting to check if need to go into more detail with checking blocks."
    if unsorted_hits == sorted_hits:
        return False
    if sorted(sorted_hits, key=lambda x: (x.ctg_pos, x.read_pos), reverse=True) == unsorted_hits:
        return False
    return True

def print_paf(outfile: io.TextIOWrapper, accepted_contigs: list, read_len: int, read_name: str,
              scaffolds: dict, k: int):
    "Print the given read mappings in PAF-like format"
    for ctg in accepted_contigs:
        ctg_run = accepted_contigs[ctg]
        sorted_mx_positions = sorted(ctg_run.hits, key=lambda x: (x.ctg_pos, x.read_pos))
        check_mapped_blocks = check_if_must_check_mapped_blocks(ctg_run.hits, sorted_mx_positions)
        if check_mapped_blocks:
            mapped_blocks = get_mapped_blocks(sorted_mx_positions)
        else:
            mapped_blocks = [sorted_mx_positions]
        for mapping in mapped_blocks:
            first_mx_mapping = mapping[0]
            last_mx_mapping = mapping[-1]
            strand_counter_list = [hit.ctg_strand == hit.read_strand for hit in mapping]
            strand_count_true = strand_counter_list.count(True)
            if strand_count_true / len(strand_counter_list) * 100 >= 50:
                strand = "+"
            else:
                strand = "-"
            target_start, target_end = min(first_mx_mapping.ctg_pos, last_mx_mapping.ctg_pos), \
                                       max(first_mx_mapping.ctg_pos, last_mx_mapping.ctg_pos) + k
            query_start, query_end = min(first_mx_mapping.read_pos, last_mx_mapping.read_pos), \
                                     max(first_mx_mapping.read_pos, last_mx_mapping.read_pos) + k
            assert query_start < query_end
            assert query_start >= 0
            assert query_end <= read_len

            out_str = f"{read_name}\t{read_len}\t{query_start}\t{query_end}\t{strand}\t" \
                      f"{ctg_run.contig}\t{scaffolds[ctg_run.contig].length}\t" \
                      f"{target_start}\t{target_end}\t{len(mapping)}\t" \
                      f"{target_end - target_start}\t255\n"
            outfile.write(out_str)
