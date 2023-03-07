from bioinformatics_textbook.inout import RosalindSolution
from bioinformatics_textbook.ch02.motif import Motif


class BA2A(RosalindSolution):
    def _solve_problem(self) -> set:
        k_d_motifs = Motif().find_k_d_motifs(
            kmer_length=self.dataset.kmer_length,
            num_allowed_mismatches=self.dataset.num_allowed_mismatches,
            dnas=self.dataset.dnas,
        )

        return self._format_rosalind_answer(k_d_motifs)
