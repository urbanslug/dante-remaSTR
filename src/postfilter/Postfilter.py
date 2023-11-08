import argparse

from ..annotation import Annotation


class PostFilter:
    """
    Class that encapsulates post-filtering.
    """

    def __init__(self, args: argparse.Namespace):
        """
        Initialize post-filter class.
        :param args: arguments of the program
        """
        self.min_flank_len = args.min_flank_len
        self.min_rep_len = args.min_rep_len
        self.min_rep_cnt = args.min_rep_cnt
        self.max_rel_error = args.max_rel_error
        self.max_abs_error = args.max_abs_error

    def quality_annotation(self, ann: Annotation, module_number: int, both_primers: bool = True) -> bool:
        """
        Is this annotation good?
        :param ann: Annotation - annotation to be evaluated
        :param module_number: int - module number
        :param both_primers: bool - do we require both primers to be present
        :return: bool - quality annotation?
        """
        is_right = ann.is_annotated_right()

        primers = ann.primers(module_number)
        has_primers = primers == 2 if both_primers else primers >= 1

        errors = (ann.has_less_errors(self.max_rel_error, relative=True) and
                  ann.has_less_errors(self.max_abs_error, relative=False))

        left_flank = sum(ann.module_bases[module_number + 1:]) >= self.min_flank_len
        right_flank = sum(ann.module_bases[:module_number]) >= self.min_flank_len
        flanks = left_flank and right_flank if both_primers else left_flank or right_flank

        repetitions = (ann.module_bases[module_number] >= self.min_rep_len and
                       ann.module_repetitions[module_number] >= self.min_rep_cnt)

        return is_right and has_primers and errors and flanks and repetitions

    def get_filtered_list(self, annotations: list[Annotation], module_number: list[int],
                          both_primers: list[bool] = None) -> tuple[list[Annotation], list[Annotation]]:
        """
        Get filtered annotations (list of modules).
        :param annotations: list(Annotation) - annotations
        :param module_number: list(int) - module numbers
        :param both_primers: list(bool) or None - do we require both primers to be present
        :return: list(Annotation), list(Annotation) - quality annotations, non-quality annotations
        """
        # adjust input if needed
        if both_primers is None:
            both_primers = [True] * len(module_number)
        assert len(both_primers) == len(module_number)

        # filter annotations
        quality_annotations = [an for an in annotations if
                               all([self.quality_annotation(an, mn, both_primers=bp) for mn, bp in zip(module_number, both_primers)])]
        filtered_annotations = [an for an in annotations if an not in quality_annotations]

        return quality_annotations, filtered_annotations

    def get_filtered(self, annotations: list[Annotation], module_number: int,
                     both_primers: bool = True) -> tuple[list[Annotation], list[Annotation]]:
        """
        Get filtered annotations.
        :param annotations: list(Annotation) - annotations
        :param module_number: int - module number
        :param both_primers: bool - do we require both primers to be present
        :return: list(Annotation), list(Annotation) - quality annotations, non-quality annotations
        """
        # pick quality annotations
        quality_annotations = [an for an in annotations if self.quality_annotation(an, module_number, both_primers=both_primers)]
        filtered_annotations = [an for an in annotations if an not in quality_annotations]

        return quality_annotations, filtered_annotations
