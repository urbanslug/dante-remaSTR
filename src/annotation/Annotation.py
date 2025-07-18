from __future__ import annotations

import re
from typing import Optional

from src.annotation import Motif

# define base mapping to regex for nucleotide symbols
base_mapping = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'R': '[GA]',
    'Y': '[CT]',
    'K': '[GT]',
    'M': '[AC]',
    'S': '[GC]',
    'W': '[AT]',
    'D': '[GAT]',
    'H': '[ACT]',
    'V': '[GCA]',
    'N': '[ACTG]'
}


class Annotation:
    """
    Encapsulate sequence of states from HMM and provide its readable representation and filters
    """

    def __init__(self, read_id: str, mate_order: int, read_seq: str, expected_seq: str, states: str, probability: float, motif: Motif):
        """
        :param read_id: str - read ID
        :param read_seq: str - read sequence
        :param mate_order: int - mate order (0 - unpaired, 1 - left pair, 2 - right pair)
        :param expected_seq: str - expected sequence as in motif
        :param states: str - sequence of states (numbers of modules)
        :param probability: Probability of generating sequence by the most likely sequence of HMM states
        :param motif: Sequence of tuples (sequence, repeats) as specified by user
        """

        # Store arguments into instance variables
        self.read_id = read_id
        self.mate_order = mate_order
        self.ordered = mate_order > 0
        self.left_pair = mate_order == 1
        self.read_seq = read_seq
        self.expected_seq = expected_seq
        self.states = states
        self.probability = probability
        self.motif = motif
        self.n_modules = len(motif.get_modules())

        # Calculate insertion/deletion/mismatch string
        self.mismatches_string = self.__get_errors()

        # Calculate number of insertions, deletions and normal bases
        self.n_insertions = self.mismatches_string.count('I')
        self.n_deletions = self.mismatches_string.count('D')
        self.n_mismatches = self.mismatches_string.count('M')

        # Number of STR motif repetitions and sequences of modules
        self.module_bases = self.__get_bases_per_module()
        self.module_repetitions = self.__get_module_repetitions()
        self.module_sequences = self.__get_module_sequences()

        # get left flank length
        self.left_flank_len = self.__get_left_flank()

    def __str__(self) -> str:
        """
        Return the annotation.
        :return: str - annotation
        """
        return '\n'.join([f'{self.read_id} {str(self.module_bases)} {str(self.module_repetitions)}', self.read_seq,
                          self.expected_seq, self.states, self.mismatches_string])

    def __get_errors(self) -> str:
        """
        Count errors in annotation and the error line.
        :return: str - error line
        """
        err_line = []
        for exp, read in zip(self.expected_seq.upper(), self.read_seq.upper()):
            if exp == '-' or read in base_mapping.get(exp, ''):
                err_line.append('_')
            elif read == '_':
                err_line.append('D')
            elif exp == '_':
                err_line.append('I')
            else:
                err_line.append('M')

        return ''.join(err_line)

    def __get_bases_per_module(self) -> tuple[int, ...]:
        """
        List of integers, each value corresponds to number of bases of input sequence that were generated by the module
        :return: Number of bases generated by each module
        """
        # Count the module states
        return tuple(self.states.count(chr(ord('0') + i)) for i in range(self.n_modules))

    def __get_left_flank(self) -> int:
        """
        Get length of a left flank.
        :return: int - number of bases of left flank before module '0' (usually module '0' is still left flank)
        """
        for i, state in enumerate(self.states):
            if state != '-':
                return i
        return len(self.states)

    def __get_module_repetitions(self) -> tuple[int, ...]:
        """
        List of integers, each value corresponds to number of repetitions of module in annotation
        :return: Number of repetitions generated by each module
        """
        # Count the module states
        repetitions = self.__get_bases_per_module()

        # Divide by the module length where applicable
        # TODO this is not right for grey ones, where only closed ones should be counted, so round is not right.
        return tuple([1 if reps == 1 and cnt > 0 else round(cnt / len(seq)) for (seq, reps), cnt in zip(self.motif.get_modules(), repetitions)])

    def __get_module_sequences(self) -> tuple[str, ...]:
        """
        List of sequences, each per module
        :return: list(str)
        """
        sequences = [''] * self.n_modules
        for i in range(self.n_modules):
            state_char = chr(ord('0') + i)
            first = self.states.find(state_char)
            if first > -1:
                last = self.states.rfind(state_char)
                sequences[i] = self.read_seq[first:last + 1]
        return tuple(sequences)

    def get_module_errors(self, module_num: int, overhang: int = None) -> tuple[int, int, int]:
        """
        Get the number of insertions and deletions or mismatches in a certain module.
        If overhang is specified, look at specified number of bases around the module as well.
        :param module_num: int - 0-based module number to count errors
        :param overhang: int - how long to look beyond module, if None, one length of STR module
        :return: int, int, int - number of insertions and deletions, mismatches, length of the interval
        """
        # get overhang as module length
        if overhang is None:
            seq, _ = self.motif.get_modules()[module_num]
            overhang = len(seq)

        # define module character
        char_to_search = chr(ord('0') + module_num)

        # if the annotation does not have this module, return 0
        if char_to_search not in self.states:
            return 0, 0, 0

        # search for the annotation of the module
        start = max(0, self.states.find(char_to_search) - overhang)
        end = min(self.states.rfind(char_to_search) + overhang + 1, len(self.states))

        # count errors
        indels = self.mismatches_string[start:end].count('I') + self.mismatches_string[start:end].count('D')
        mismatches = self.mismatches_string[start:end].count('M')

        # return indels, mismatches, and length
        return indels, mismatches, end - start

    def info_value(self) -> (int, int):
        """
        Evaluate the info value of the Annotation.
        :return: int, int
        """
        return sum(self.module_repetitions), sum(self.module_bases)

    def info_value_str(self, index_str: int) -> (int, int, int):
        """
        Evaluate the info value of the Annotation at the STR location
        :param index_str: int - index of the STR
        :return: int, int, int
        """
        return self.primers(index_str), self.module_repetitions[index_str], self.module_bases[index_str]

    def has_required_modules(self, required_repetitions: list[int]) -> bool:
        """
        Validate, if read is annotated with sufficient number of modules
        :param required_repetitions: list of number of required repetitions, one for each module
        :return: True, if annotation has required number of repetition for each module
        """
        if not required_repetitions:
            return True

        for repetition, required_repetition in zip(self.module_repetitions, required_repetitions):
            if repetition < required_repetition:
                return False
        return True

    def has_required_bases(self, required_bases: list[int]) -> bool:
        """
        Validate, if read bases are sufficiently annotated
        :param required_bases: list of number of required annotated bases, one for each module
        :return: True, if annotation has required number of annotated bases for each module
        """
        if not required_bases:
            return True

        for bases, required_base in zip(self.module_bases, required_bases):
            if bases < required_base:
                return False
        return True

    def has_one_primer(self, required_bases: list[int], required_repetitions: list[int], index_rep: int, index_rep2: int = None) -> bool:
        """
        Validate, if at least one primer is sufficiently annotated, in case of 2 repetitions, we check only the
        :param required_bases: list of number of required annotated bases, one for each module
        :param required_repetitions: list of number of required annotated modules
        :param index_rep: int - index of the repetition, that we are looking at
        :param index_rep2: int - index of the second repetition, that we are looking at
        :return: True, if annotation has at least one primer is sufficiently annotated
        """
        # if it is not interesting just throw it away
        if not self.is_annotated_right():
            return False

        # if no filter, accept all
        if required_bases is None and required_repetitions is None:
            return True

        def check_reqs(index: int) -> bool:
            """
            Check if requirements are met on 'index'
            :param index: int - index of a module
            :return: bool - requirements are met?
            """
            if required_bases is not None and self.module_bases[index] < required_bases[index]:
                return False
            if required_repetitions is not None and self.module_repetitions[index] < required_repetitions[index]:
                return False
            return True

        if index_rep2 is None:
            index_rep2 = index_rep
        index_rep, index_rep2 = min(index_rep, index_rep2), max(index_rep, index_rep2)

        # if the module has right primer, and it is clipped on the left
        if index_rep2 + 1 < len(required_bases) and self.states[0] != '-' and check_reqs(index_rep2 + 1) and check_reqs(index_rep):
            return True
        # if the module has left primer, and it is clipped on the right
        if index_rep - 1 >= 0 and self.states[-1] != '-' and check_reqs(index_rep - 1) and check_reqs(index_rep2):
            return True
        return False

    def has_less_errors(self, max_errors: float | int, relative=False) -> bool:
        """
        Check if this annotation has fewer errors than max_errors. Make it relative to the annotated length if relative is set.
        :param max_errors: int/float - number of max_errors (relative if relative is set)
        :param relative: bool - if the errors are relative to the annotated length
        :return: bool - True if the number of errors is less than allowed
        """
        errors = self.n_deletions + self.n_insertions + self.n_mismatches

        if max_errors is None or errors == 0:
            return True

        if relative:
            return errors / float(sum(self.module_bases)) <= max_errors
        else:
            return errors <= max_errors

    def primers(self, index_rep: int) -> int:
        """
        Count how any primers it has on repetition index.
        :param index_rep: int - index of the repetition, that we are looking at
        :return: int - number of primers (0-2)
        """
        primers = 0
        if index_rep > 0 and self.module_repetitions[index_rep - 1] > 0:
            primers += 1
        if index_rep + 1 < len(self.module_repetitions) and self.module_repetitions[index_rep + 1] > 0:
            primers += 1
        return primers

    def is_annotated_right(self) -> bool:
        """
        Is it annotated in a way that it is interesting? More than one module annotated + modules are not missing in the middle.
        :return: bool - annotated right?
        """

        # remove those that starts/ends in background but don't have a neighbour module
        starts_background = self.states[0] in '_-'
        ends_background = self.states[-1] in '_-'
        if starts_background and self.module_repetitions[0] == 0:
            return False
        if ends_background and self.module_repetitions[-1] == 0:
            return False

        # remove those with jumping modules
        started = False
        ended = False
        for repetition in self.module_repetitions:
            if repetition > 0:
                started = True
                if ended:
                    return False
            if repetition == 0 and started:
                ended = True

        # pass?
        return True

    def same_start_fragment(self, annotation: Annotation) -> bool:
        """
        Return True if both sequences can be produced from the same start of a fragment.
        :param annotation: Annotation - second annotation
        :return: bool - True if both sequences can be produced from the same start of a fragment
        """
        comp_length = min(len(self.read_seq), len(annotation.read_seq))
        return self.read_seq[:comp_length] == annotation.read_seq[:comp_length]

    def same_end_fragment(self, annotation: Annotation) -> bool:
        """
        Return True if both sequences can be produced from the same end of a fragment.
        :param annotation: Annotation | None - second annotation
        :return: bool - True if both sequences can be produced from the same end of a fragment
        """
        comp_length = min(len(self.read_seq), len(annotation.read_seq))
        return self.read_seq[-comp_length:] == annotation.read_seq[-comp_length:]

    def get_str_repetitions(self, index_str: int) -> Optional[tuple[bool, int]]:
        """
        Get the number of str repetitions for a particular index.
        :param index_str: int - index of a str
        :return: (bool, int) - closed?, number of str repetitions
        """
        if self.is_annotated_right():
            primer1 = index_str > 0 and self.module_repetitions[index_str - 1] > 0
            primer2 = index_str + 1 < len(self.module_repetitions) and self.module_repetitions[index_str + 1] > 0
            if primer1 or primer2:
                return primer1 and primer2, self.module_repetitions[index_str]
        return None

    @staticmethod
    def find_with_regex(read_sequence: str, motif_sequence: str, search_pos: int = 0) -> int:
        """
        Find the first occurrence of a motif sequence in the read sequence using regular expressions.
        :param read_sequence: The sequence to search in.
        :param motif_sequence: The motif sequence (as a regex) to search for.
        :param search_pos: The position to start the search from.
        :return: int - The start position of the first occurrence of the motif sequence. Returns -1 if not found.
        """
        # convert motif sequence to regex
        motif_regex = ''.join(base_mapping[char] for char in motif_sequence)

        # compile the regular expression pattern
        pattern = re.compile(motif_regex)

        # search for the pattern in the read sequence starting from search_pos
        match = pattern.search(read_sequence, search_pos)

        # return the start position if a match is found, else return -1
        return match.start() if match else -1

    def get_nomenclature(self, index_rep: int = None, index_rep2: int = None, include_flanking: bool = True) -> str:
        """
        Get HGVS nomenclature.
        :param index_rep: int - index of the first repetition (None if include all)
        :param index_rep2: int - index of the second repetition (None if include all)
        :param include_flanking: boolean - include flanking regions (i.e. first and last module)
        :return: str - HGVS nomenclature string
        """
        # prepare data
        if index_rep is not None:
            if index_rep2 is not None:
                data = zip([self.module_repetitions[index_rep], self.module_repetitions[index_rep2]], [self.motif[index_rep], self.motif[index_rep2]],
                           [self.module_sequences[index_rep], self.module_sequences[index_rep2]])
            else:
                data = zip([self.module_repetitions[index_rep]], [self.motif[index_rep]], [self.module_sequences[index_rep]])
        elif include_flanking:
            data = zip(self.module_repetitions, self.motif.get_modules(), self.module_sequences)
        else:
            data = zip(self.module_repetitions[1:-1], self.motif[1:-1], self.module_sequences[1:-1])

        # iterate and build the nomenclature string
        nomenclatures = []
        for repetitions, (motif_sequence, _), read_sequence in data:
            nomenclature = ''
            if repetitions == 1:
                if len(read_sequence) > 0:
                    nomenclature += f'{read_sequence}[1]'
            else:
                reps = 0
                search_pos = 0
                found_rep_seq = ''
                while True:
                    search_found = self.find_with_regex(read_sequence, motif_sequence, search_pos)
                    if search_found == search_pos:
                        # setup current rep. sequence
                        if reps == 0:
                            found_rep_seq = read_sequence[search_found:search_found + len(motif_sequence)]

                        if read_sequence[search_found:search_found + len(motif_sequence)] == found_rep_seq:
                            # regular continuation
                            reps += 1
                        else:
                            # interruption, but in line with searched motif
                            nomenclature += f'{found_rep_seq}[{reps}]'
                            found_rep_seq = read_sequence[search_found:search_found + len(motif_sequence)]
                            reps = 1
                    elif search_found == -1:
                        # the end, we did not find any other STRs
                        if reps > 0:
                            nomenclature += f'{found_rep_seq}[{reps}]'
                        if len(read_sequence[search_pos:]) > 0:
                            nomenclature += f'{read_sequence[search_pos:]}[1]'
                        break
                    else:
                        # some interruption
                        if reps > 0:
                            nomenclature += f'{found_rep_seq}[{reps}]'
                        if len(read_sequence[search_pos:search_found]) > 0:
                            nomenclature += f'{read_sequence[search_pos:search_found]}[1]'
                        found_rep_seq = read_sequence[search_found:search_found + len(motif_sequence)]
                        reps = 1
                    # update search pos and iterate
                    search_pos = search_found + len(motif_sequence)
            nomenclatures.append(nomenclature)

        return '\t'.join(nomenclatures)

    def get_shortened_annotation(self, shorten_length: int) -> Annotation:
        """
        Get shortened annotation with specified shorten length beyond annotated modules.
        :param shorten_length: int - how many bases to keep beyond modules
        :return: Annotation - shortened annotation
        """

        # search for start
        start = -1
        for i in range(len(self.states)):
            if self.states[i] != '-':
                start = i
                break

        # search for end
        end = -1
        for i in range(len(self.states) - 1, -1, -1):
            if self.states[i] != '-':
                end = i
                break

        # adjust start and end for shorten length
        start = max(start - shorten_length, 0)
        end = min(end + 1 + shorten_length, len(self.states))  # +1 for use as list range

        # return shortened Annotation
        return Annotation(self.read_id, self.mate_order, self.read_seq[start:end], self.expected_seq[start:end],
                          self.states[start:end], self.probability, self.motif)
