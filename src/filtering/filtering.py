import numpy as np


def has_good_quality(row: dict, min_qual: int, first_module: int, last_module: int) -> bool:
    """
    Checks if the read has good quality everywhere in the annotated part.
    :param row: dict - the read, its sequence, data, quality?
    :param min_qual: int - minimal quality to have
    :param first_module: int - first module of the annotation to look at
    :param last_module: int - last module of the annotation to look at
    :return: bool - True if the read has enough quality bases in the annotated part
    """
    # first adjust the lengths due to insertions
    indices_to_except = set([i for i, ltr in enumerate(row['read']) if ltr == '_'])
    seq = ''.join(c for i, c in enumerate(row['read']) if i not in indices_to_except)
    ref = ''.join(c for i, c in enumerate(row['reference']) if i not in indices_to_except)
    mod = ''.join(c for i, c in enumerate(row['modules']) if i not in indices_to_except)
    qual = [ord(q) - ord('!') for q in row['quality']]

    assert len(seq) == len(ref) == len(mod) == len(qual), (row['motif'], row['read_sn'], len(seq), len(qual), row['read'], seq, row['quality'])

    # identify annotated place
    annot_start = mod.find(chr(ord('0') + first_module))
    if annot_start == -1:
        if mod[0] != '-' and ord(mod[0]) - ord('0') <= last_module:
            annot_start = 0
        else:
            return True
    annot_end = mod.rfind(chr(ord('0') + last_module))
    annot_end = len(mod) if annot_end == -1 else annot_end + 1

    # look if there is enough quality in the annotated part
    return all(q >= min_qual for q in qual[annot_start:annot_end])


def cut_low_quality(row: dict, min_qual: int) -> dict:
    """
    Cut parts of the read with low quality
    :param row: dict - the read, its sequence, data, quality?
    :param min_qual: int - minimal quality to have
    :return: dict - row with updated sequences/qualities
    """
    # first adjust the quality column due to insertions
    indices_ins = [i for i, ltr in enumerate(row['read']) if ltr == '_']
    result = []
    last_idx = 0
    for idx in sorted(indices_ins):
        result.append(row['quality'][last_idx:idx])  # Append substring up to current index
        result.append('X')  # Insert 'X' (high quality)
        last_idx = idx  # Update last index to current
    result.append(row['quality'][last_idx:])  # Append the remaining part of the string
    qual = [ord(q) - ord('!') for q in ''.join(result)]

    assert len(row['read']) == len(qual), (row['motif'], row['read_sn'], len(qual), row['read'], row['quality'])

    # cut low quality from relevant parts of row
    keep_ids = np.array(qual) >= min_qual
    row['read'] = ''.join(np.array(list(row['read']))[keep_ids])
    row['reference'] = ''.join(np.array(list(row['reference']))[keep_ids])
    row['modules'] = ''.join(np.array(list(row['modules']))[keep_ids])
    row['quality'] = ''.join([chr(q + ord('!')) for q in np.array(qual)[keep_ids]])

    # return the modified row
    return row
