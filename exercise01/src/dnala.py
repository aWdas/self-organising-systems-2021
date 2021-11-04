import math

from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.SeqRecord import SeqRecord

from .util import dist, generalize, seq_str


def distance_matrix(sequences: [SeqRecord]) -> [[int]]:
    return [[dist(seq_str(seq_a), seq_str(seq_b)) for seq_b in sequences] for seq_a in sequences]


def dnala(sequences: [SeqRecord]) -> [str]:
    print("Calculating distance matrix")
    dist_matrix = distance_matrix(sequences)

    generalized_sequences = []
    generalized_distance = 0
    remaining_sequences = list(range(len(sequences)))

    def closest_sequences(i):
        min_dist = min(dist_matrix[i][j] if i != j else math.inf for j in remaining_sequences)
        return [j for j in remaining_sequences if dist_matrix[i][j] == min_dist and j != i]

    while len(remaining_sequences) > 0:
        if (len(remaining_sequences)) == 3:
            [x, y, z] = remaining_sequences
            intermediate_sequence = generalize(seq_str(sequences[x]), seq_str(sequences[y]))
            generalized_sequences.append(generalize(intermediate_sequence, seq_str(sequences[z])))
            generalized_distance = generalized_distance + \
                                   dist_matrix[x][y] + \
                                   dist(intermediate_sequence, seq_str(sequences[z]))
            remaining_sequences = []

        for s in remaining_sequences:
            candidates = closest_sequences(s)
            match_found = False
            for c in candidates:
                reverse_candidates = closest_sequences(c)
                if s in reverse_candidates:
                    generalized_sequence = generalize(seq_str(sequences[s]), seq_str(sequences[c]))
                    print(f"Matched {sequences[s]} to {sequences[c]}: {generalized_sequence}")
                    generalized_sequences.append(generalized_sequence)
                    remaining_sequences.remove(s)
                    remaining_sequences.remove(c)
                    generalized_distance = generalized_distance + dist_matrix[s][c]
                    match_found = True
                    break
            if match_found:
                break

    return generalized_sequences, generalized_distance
