from bioinformatics_textbook.inout import RosalindSolution
from bioinformatics_textbook.ch02.median_string import MedianString
from bioinformatics_textbook.ch02.motif import Motif


class BA2A(RosalindSolution):
    def _solve_problem(self) -> set:
        k_d_motifs = Motif().find_k_d_motifs(
            kmer_length=self.dataset.kmer_length,
            num_allowed_mismatches=self.dataset.num_allowed_mismatches,
            dnas=self.dataset.dnas,
        )

        return self._format_rosalind_answer(k_d_motifs)


class BA2B(RosalindSolution):
    def _solve_problem(self) -> int:
        median_strings = MedianString().find_median_strings(
            kmer_length=self.dataset.k,
            dnas=self.dataset.dnas,
        )

        return median_strings[0]


class BA2H(RosalindSolution):
    def _solve_problem(self) -> int:
        distance = MedianString().compute_pattern_strings_distance(
            pattern=self.dataset.pattern,
            dnas=self.dataset.dnas,
        )

        return distance
