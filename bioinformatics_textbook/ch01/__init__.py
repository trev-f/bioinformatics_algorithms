from bioinformatics_textbook.inout import (
    RosalindSolution
)

from bioinformatics_textbook.ch01.pattern_occurrences import (
    PatternOccurrences
)

from bioinformatics_textbook.ch01.frequent_words import (
    FrequentWords
)

from bioinformatics_textbook.ch01.reverse_complement import ReverseComplement


class BA1A(RosalindSolution):
    def _solve_problem(self) -> str:
        kmer_count = PatternOccurrences().count_pattern(
            text=self.dataset.text,
            pattern=self.dataset.pattern
        )

        return kmer_count


class BA1B(RosalindSolution):
    def _solve_problem(self) -> str:
        most_freq_words = FrequentWords().find_most_freq_words(
            text=self.dataset.text,
            kmer_length=self.dataset.kmer_length,
        )

        return self._format_rosalind_answer(most_freq_words)


class BA1C(RosalindSolution):
    def _solve_problem(self) -> str:
        rev_comp = ReverseComplement().reverse_complement_dna(dna=self.dataset.pattern)

        return rev_comp


class BA1D(RosalindSolution):
    def _solve_problem(self) -> str:
        starting_positions = PatternOccurrences().find_starting_positions(pattern=self.dataset.pattern, genome=self.dataset.genome)
        starting_positions = self._convert_iterable_to_list_of_str(starting_positions)
        
        return self._format_rosalind_answer(starting_positions)


class BA1I(RosalindSolution):
    def _solve_problem(self) -> str:
        most_freq_words = FrequentWords().find_most_freq_words_with_mismatches(
            text=self.dataset.text,
            kmer_length=self.dataset.kmer_length,
            num_allowed_mismatches=self.dataset.hamming_dist,
        )

        return self._format_rosalind_answer(most_freq_words)


class BA1J(RosalindSolution):
    def _solve_problem(self) -> str:
        most_freq_words = FrequentWords().find_most_freq_words_with_mismatches_and_rc(
            text=self.dataset.text,
            kmer_length=self.dataset.kmer_length,
            num_allowed_mismatches=self.dataset.hamming_dist,
        )

        return self._format_rosalind_answer(most_freq_words)


class BA1N(RosalindSolution):
    def _solve_problem(self) -> str:
        neighborhood = FrequentWords().find_neighbors(
            pattern=self.dataset.pattern,
            num_allowed_mismatches=self.dataset.hamming_dist,
        )

        return self._format_rosalind_answer(neighborhood, sep='\n')
