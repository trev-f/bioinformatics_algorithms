from bioinformatics_textbook.dna import DNA

class MedianString:
    def find_median_string(self) -> None:
        pass

    def compute_pattern_strings_distance(self, pattern: DNA, dnas: list[DNA]) -> int:
        """Compute the distance between a pattern and a collection of DNA sequences

        :param pattern: k-mer
        :type pattern: DNA
        :param dnas: Collection of DNA sequences
        :type dnas: list[DNA]
        :return: Distance between pattern and DNA sequences
        :rtype: int
        """
        global_distance = 0
        for dna in dnas:
            local_distance = len(pattern)
            for kmer in dna.generate_kmers(kmer_length=len(pattern)):
                hamming_distance = pattern.compute_hamming_distance(dna_q=kmer)
                local_distance = hamming_distance if hamming_distance < local_distance else local_distance
            
            global_distance += local_distance

        return global_distance
