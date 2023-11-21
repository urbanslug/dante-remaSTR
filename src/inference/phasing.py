import os
from collections import Counter

from src.annotation import Annotation


def phase(annotations: list[Annotation], module_number1: int, module_number2: int) -> tuple[tuple[str, str], tuple[str, str, str]]:
    """
    Infer phasing based on the Annotations.
    :param annotations: list(Annotation) - good (blue) annotations
    :param module_number1: int - index of a repetition
    :param module_number2: int - index of the second repetition
    :return: tuple - predicted symbols and confidences
    """
    # resolve trivial case
    if len(annotations) == 0:
        return ('-/-', '-/-'), ('-/0', '-/0', '-/0')

    # gather module repetitions from annotations and count them
    repetitions = Counter([(ann.module_repetitions[module_number1], ann.module_repetitions[module_number2]) for ann in annotations])

    # pick the highest two
    most_common = repetitions.most_common(2)
    rep1, cnt1 = most_common[0]
    rep2, cnt2 = most_common[1] if len(most_common) >= 2 else (('-', '-'), 0)

    # output phasing with number of supported reads
    phasing = (f'{rep1[0]}/{rep1[1]}', f'{rep2[0]}/{rep2[1]}')
    supported_reads = (f'{cnt1 + cnt2}/{len(annotations)}', f'{cnt1}/{len(annotations)}', f'{cnt2}/{len(annotations)}')
    return phasing, supported_reads


def save_phasing(phasing_file: str, phasing: tuple[str, str], supp_reads: tuple[str, str, str]) -> None:
    """
    Save the phasing information to a file.
    :param phasing_file: str - filename where to save the phasing information
    :param phasing: tuple - phasing information
    :param supp_reads: tuple - supporting reads (number of reads that support the
    """
    # save into a file:
    with open(phasing_file, 'w') as f:
        f.write(f'{phasing[0]}\t{phasing[1]}\n')
        f.write(f'{supp_reads[0]}\t{supp_reads[1]}\t{supp_reads[2]}\n')


def load_phasing(phasing_file: str) -> tuple[tuple[str, str], tuple[str, str, str]] | None:
    """
    Load the phasing information from a file.
    :param phasing_file: str - filename from where to load the phasing information
    :return: tuple - predicted symbols and confidences
    """
    if not os.path.exists(phasing_file):
        return None

    # load the file
    with open(phasing_file, 'r') as f:
        phasing = f.readline().strip().split('\t')
        supp_reads = f.readline().strip().split('\t')

    return (phasing[0], phasing[1]), (supp_reads[0], supp_reads[1], supp_reads[2])
