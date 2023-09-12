import re
from typing import List, Tuple


class Motif:
    """
    Class to represent DNA motifs.

    :ivar chrom: Chromosome name.
    :ivar start: Start position of the motif.
    :ivar end: End position of the motif.
    :ivar modules: A list of tuples containing sequence and repetition count.
    """

    def __init__(self, motif: str):
        """
        Initialize a Motif object.
        :param motif: The motif string in the format "chrom:start_end[A][B]..."
        """
        # remove whitespace
        self.name = motif.strip().replace(' ', '')

        # extract prefix, first number, second number
        self.chrom, start, end, remainder = re.match(r'([^:]+):g\.(\d+)_(\d+)(.*)', self.name).groups()

        # extract sequence and repetition count
        self.modules = [(str(seq), int(num)) for seq, num in re.findall(r'([A-Z]+)\[(\d+)\]', remainder)]
        self.modules = [('left_flank', 1)] + self.modules + [('right_flank', 1)]

        # convert to ints
        self.start = int(start)
        self.end = int(end)


    def __getitem__(self, index: int) -> Tuple[str, int]:
        """
        Returns module at given index.
        :param index: The index of the module to fetch.
        :return: The module at the given index.
        """
        return self.modules[index]


    def __str__(self) -> str:
        """
        Returns string representation of the Motif object.
        :return: String representation in the format "chrom:start_end[A][B]..."
        """
        return f'{self.chrom}:g.{self.start}_{self.end}' + self.modules_str(include_flanks=False)


    def modules_str(self, include_flanks: bool = False) -> str:
        """
        Returns string representation of modules
        :param include_flanks: bool - include flank modules?
        :return: String representation of modules
        """
        if include_flanks:
            return ''.join([f'{seq}[{num}]' for seq, num in self.modules])
        return ''.join([f'{seq}[{num}]' for seq, num in self.modules[1:-1]])


    def dir_name(self) -> str:
        """
        Returns possible directory name of the motif.
        :return: str - directory name for the motif
        """
        return self.name.replace(':', '-').replace('.', '_')


    def get_modules(self) -> List[Tuple[str, int]]:
        """
        Returns list of modules.
        :return: List of tuples containing sequence and repetition count.
        """
        return self.modules


    def get_repeating_modules(self) -> List[Tuple[int, str, int]]:
        """
        Returns list of modules with more than one repetition.
        :return: List of tuples containing index, sequence, and repetition count.
        """
        return [(i, seq, num) for i, (seq, num) in enumerate(self.modules) if num > 1]


    def get_location_subpart(self, index: int) -> Tuple[int, int]:
        """
        Returns the chromosome location of a subpart of a motive
        :param index: int - index of a module
        :return: start and end location of the subpart
        """
        start = self.start
        for module in self.modules[1: index]:
            seq, rep = module
            start += len(seq) * rep

        return start, start + len(self.modules[index][0]) * self.modules[index][1]
