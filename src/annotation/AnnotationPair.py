from __future__ import annotations

from src.annotation import Annotation


class AnnotationPair:
    """
    Encapsulate annotation pairs.
    """

    def __init__(self, ann1: Annotation, ann2: Annotation | None):
        """
        Initialize the AnnotationPair object
        :param ann1: Annotation - first annotation of a pair
        :type ann1: Annotation
        :param ann2: Annotation - second annotation of a pair
        :type ann2: Annotation | None
        """
        assert ann1 is not None
        self.ann1 = ann1
        self.ann2 = ann2

    def __eq__(self, second_pair: AnnotationPair) -> bool:
        """
        Annotation pairs equal when they are produced by the same fragment.
        :param second_pair: AnnotationPair - second annotation pair
        :type second_pair: AnnotationPair
        :return: bool - whether the annotation pairs are produced by the same fragment
        """
        # first check if we deal with simple annotations
        if self.ann1 is None:
            return second_pair.ann2 is not None and self.ann2.same_end_fragment(second_pair.ann2)

        if self.ann2 is None:
            return second_pair.ann1 is not None and self.ann1.same_start_fragment(second_pair.ann1)

        # return the full comparison
        return ((second_pair.ann1 is None or self.ann1.same_start_fragment(second_pair.ann1)) and
                (second_pair.ann2 is None or self.ann2.same_end_fragment(second_pair.ann2)))

    def __str__(self):
        return '\n'.join(['1', str(self.ann1), '2', str(self.ann2)])

    def has_required_modules(self, required_repetitions: list[int]) -> bool:
        """
        Validate, if read modules are sufficiently annotated (in at least one of the reads)
        :param required_repetitions: list of number of required annotated modules
        :return: True, if one of the annotations has required number of annotated modules
        """
        left = self.ann1 is not None and self.ann1.has_required_modules(required_repetitions)
        right = self.ann2 is not None and self.ann2.has_required_modules(required_repetitions)
        return left or right

    def has_required_bases(self, required_bases: list[int]) -> bool:
        """
        Validate, if read bases are sufficiently annotated (in at least one of the reads)
        :param required_bases: list of number of required annotated bases, one for each module
        :return: True, if one of the annotations has required number of annotated bases for each module
        """
        left = self.ann1 is not None and self.ann1.has_required_bases(required_bases)
        right = self.ann2 is not None and self.ann2.has_required_bases(required_bases)
        return left or right

    def has_one_primer(self, required_bases: list[int], required_repetitions: list[int], index_rep: int, index_rep2: int = None) -> bool:
        """
        Validate, if at least one primer is sufficiently annotated (in at least one of the reads)
        :param required_bases: list of number of required annotated bases, one for each module
        :param required_repetitions: list of number of required annotated modules
        :param index_rep: int - index of the repetition, that we are looking at
        :param index_rep2: int - index of the second repetition, that we are looking at
        :return: True, if one of the annotations has at least one primer is sufficiently annotated
        """
        left = self.ann1 is not None and self.ann1.has_one_primer(required_bases, required_repetitions, index_rep, index_rep2)
        rght = self.ann2 is not None and self.ann2.has_one_primer(required_bases, required_repetitions, index_rep, index_rep2)
        return left or rght

    def more_info_than(self, second_pair: AnnotationPair) -> bool:
        """
        Check if this AnnotationPair has more info than second AnnotationPair
        :param second_pair: AnnotationPair - second annotation pair
        :return: bool - True if it has more info than the other AnnotationPair
        """

        # first compare if both has annotations:
        if self.ann2 is None:
            return False

        if second_pair.ann2 is None:
            return True

        # then return those that have more annotated modules, or bases:
        m1f, b1f = self.ann1.info_value()
        m2f, b2f = self.ann2.info_value()

        m1s, b1s = second_pair.ann1.info_value()
        m2s, b2s = second_pair.ann2.info_value()

        return (m1f + m2f, b1f + b2f) > (m1s + m2s, b1s + b2s)

    def get_more_informative_annotation(self, index_str: int = None) -> Annotation:
        """
        Return the more informative Annotation from the two.
        :return: Annotation - the more informative annotation
        """
        # If the second is non-existent, return the first
        if self.ann2 is None:
            return self.ann1

        # It the STR index is not provided return the "absolutely" more informative, else return the more informative on that STR
        if index_str is None:
            if self.ann1.info_value() >= self.ann2.info_value():
                return self.ann1
            else:
                return self.ann2
        else:
            if self.ann1.info_value_str(index_str) >= self.ann2.info_value_str(index_str):
                return self.ann1
            else:
                return self.ann2

    def get_str_repetitions(self, index_str: int) -> (bool, int):
        """
        Get the number of str repetitions for a particular index.
        :param index_str: int - index of a str
        :return: (bool, int) - closed?, number of str repetitions
        """

        def add_primers(annotation: Annotation, module_num: int) -> list[(bool, int)]:
            """
            Add str repetitions from one annotation into an array of results.
            :param annotation: Annotation - input annotation
            :param module_num: int - index of a str
            :return: list(bool, int) - closed?, number of str repetitions
            """
            primer1 = module_num > 0 and annotation.module_repetitions[module_num - 1] > 0
            primer2 = module_num + 1 < len(annotation.module_repetitions) and annotation.module_repetitions[module_num + 1] > 0
            if primer1 or primer2:
                return [(primer1 and primer2, annotation.module_repetitions[module_num])]
            return []

        results = []
        if self.ann1 is not None and self.ann1.is_annotated_right():
            results.extend(add_primers(self.ann1, index_str))
        if self.ann2 is not None and self.ann2.is_annotated_right():
            results.extend(add_primers(self.ann2, index_str))

        # return the highest (first see if it is closed and then pick the highest number):
        if len(results) > 0:
            return sorted(results)[-1]
        return None


def annotations_to_pairs(annots: list[Annotation]) -> list[AnnotationPair]:
    """
    Convert an array of annotations to annotation pairs array.
    :param annots: list(Annotation) - annotations
    :return: list(AnnotationPair)
    """

    # sort:
    sorted_list = sorted(annots, key=lambda annot: (annot.read_id, annot.left_pair))

    # remove duplicates:
    seen = set()
    deduplicated = []
    for ann in sorted_list:
        if (ann.read_id, ann.left_pair) not in seen:
            seen.add((ann.read_id, ann.left_pair))
            deduplicated.append(ann)

    result = []

    i = 0
    while i < len(deduplicated):
        if i + 1 < len(deduplicated) and deduplicated[i].read_id == deduplicated[i + 1].read_id:
            assert deduplicated[i].ordered and deduplicated[i + 1].ordered and deduplicated[i].left_pair != deduplicated[i + 1].left_pair
            result.append(AnnotationPair(deduplicated[i], deduplicated[i + 1]))
            i += 2
        else:
            result.append(AnnotationPair(deduplicated[i], None))
            i += 1

    return result


def pairs_to_annotations(annotation_pairs: list[AnnotationPair]) -> list[Annotation]:
    """
    Convert an array of annotations pairs to annotation array.
    :param annotation_pairs: list(AnnotationPair) - annotations
    :return: list(Annotation)
    """
    annots = []
    for ap in annotation_pairs:
        if ap.ann1 is not None:
            annots.append(ap.ann1)
        if ap.ann2 is not None:
            annots.append(ap.ann2)

    return annots


def pairs_to_annotations_pick(annotation_pairs: list[AnnotationPair], index_str: int | None) -> list[Annotation]:
    """
    Convert an array of annotations pairs to annotation array. Leave only the more informative one.
    :param annotation_pairs: list(AnnotationPair) - annotations
    :param index_str: int: index of the STR to look at or None if we should get the whole motif into account
    :return: list(Annotation)
    """
    return [ap.get_more_informative_annotation(index_str) for ap in annotation_pairs]


# TODO make this faster/parallelize - takes too long when number of found reads is more than 10.000-100.000 (1min at 7.500, 4h at 400.000)
def remove_pcr_duplicates(annot_pairs: list[AnnotationPair]) -> (list[AnnotationPair], list[AnnotationPair]):
    """
    Remove PCR duplicates -- deduplicate the annotation pair list.
    :param annot_pairs: list(AnnotationPair) - list of Annotation Pairs
    :return: list(AnnotationPair), list(AnnotationPair) - deduplicated list and duplications
    """
    def remove_none(ann_pairs: list[AnnotationPair], first: bool = True) -> list[AnnotationPair]:
        """
        Remove annotation pairs with None first (second) annotation from the pair.
        :param ann_pairs: list(AnnotationPair) - list of Annotation Pairs
        :param first: bool - look at the first annotation from the pair?
        :return: list(AnnotationPair) - list of Annotation Pairs without None pairs
        """
        return [ap for ap in ann_pairs if (first and ap.ann1 is not None) or (not first and ap.ann2 is not None)]

    def restore_none(pairs_with_none: list[AnnotationPair], first: bool = True) -> list[AnnotationPair]:
        """
        Restore annotation pairs with None first (second) annotation from the pair.
        :param pairs_with_none: list(AnnotationPair) - list of Annotation Pairs with None annotations
        :param first: bool - look at the first annotation from the pair?
        :return: list(AnnotationPair) - list of Annotation Pairs with restored Annotation Pairs with None annotation
        """
        return [ap for ap in pairs_with_none if (first and ap.ann1 is None) or (not first and ap.ann2 is None)]

    def deduplicate(ann_pairs: list[AnnotationPair]) -> (list[AnnotationPair], list[AnnotationPair]):
        """
        Remove PCR duplicates -- deduplicate the annotation pair list.
        :param ann_pairs: list(AnnotationPair) - list of Annotation Pairs sorted by annotation_1 or annotation_2
        :return: list(AnnotationPair), list(AnnotationPair) - deduplicated list and duplications
        """
        dedup = []
        duplic = []

        if not ann_pairs:
            return [], []

        # Find duplicates by comparing neighbours in sorted list
        prev_ap = ann_pairs[0]
        for curr_ap in ann_pairs[1:]:
            if prev_ap == curr_ap:
                if prev_ap.more_info_than(curr_ap):
                    duplic.append(curr_ap)
                else:
                    duplic.append(prev_ap)
                    prev_ap = curr_ap
            else:
                dedup.append(prev_ap)
                prev_ap = curr_ap

        dedup.append(prev_ap)

        return dedup, duplic

    if not annot_pairs:
        return [], []

    # Deduplication according to first annotation in pair
    curr_pairs = remove_none(annot_pairs, True)
    curr_pairs = sorted(curr_pairs, key=lambda ap: ap.ann1.read_seq)
    deduplicated_1, duplicates_1 = deduplicate(curr_pairs)

    curr_pairs = deduplicated_1 + restore_none(annot_pairs, True)

    # Deduplication according to second annotation in pair
    curr_pairs = remove_none(curr_pairs, False)
    curr_pairs = sorted(curr_pairs, key=lambda ann: ann.ann2.read_seq[::-1])
    deduplicated_2, duplicates_2 = deduplicate(curr_pairs)

    deduplicated = deduplicated_2 + restore_none(deduplicated_1, False)
    duplicates = duplicates_1 + duplicates_2

    return deduplicated, duplicates
