import logging


logger = logging.getLogger(__name__)

class FrequentWords:
    
    def find_most_freq_words_with_mismatches(self, text: str, kmer_length: int, num_allowed_mismatches: int) -> list:
        # initialize empty collections and 
        patterns = []
        freq_table = {}
        n = len(text)
        for i in range(n - kmer_length + 1):
            pattern = text[i: i+kmer_length]
            neighborhood = self.find_neighbors(pattern, num_allowed_mismatches)
            for neighbor in neighborhood:
                # if a k-mer is not present in frequency table, add it and assign a value of 1,
                # otherwise, increment the count
                freq_table[neighbor] = freq_table.get(neighbor, 0) + 1
        
        max_freq = max(freq_table.values())
        for kmer in freq_table.keys():
            if freq_table[kmer] == max_freq:
                patterns.append(kmer)

        return patterns


    def find_neighbors(self, pattern: str, num_allowed_mismatches: int) -> list:
        neighbors = []
        
        return neighbors
