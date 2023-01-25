from bioinformatics_textbook.code_challenges.frequent_words import (
    FrequentWords
)

from bioinformatics_textbook.code_challenges.inout import (
    RosalindSolution
)


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
