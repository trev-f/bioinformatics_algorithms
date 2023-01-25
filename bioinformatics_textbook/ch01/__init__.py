from bioinformatics_textbook.ch01.pattern_count import (
    PatternCount
)

from bioinformatics_textbook.ch01.frequent_words import (
    FrequentWords
)

from bioinformatics_textbook.ch01.inout import (
    RosalindSolution
)


class BA1A(RosalindSolution):
    def _solve_problem(self) -> str:
        kmer_count = PatternCount().count_pattern(
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
