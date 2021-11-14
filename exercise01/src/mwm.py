from Bio import Align, SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
import networkx as nx
from Bio.pairwise2 import format_alignment
from networkx.algorithms.matching import max_weight_matching

from .util import dist, seq_str, generalize, build_pairwise_alignments, fix_odd_sequence_count


def mwm_exact(sequences: [SeqRecord], alignment_result):
    g = nx.Graph()
    m = 10000000000
    generalized_sequences = []
    generalized_distance = 0

    # build a map of the pairwise distances between the sequences,
    # and keep track of the pair with the minimum distance
    if alignment_result is None:
        print("No pre-computed pairwise alignments were provided -> calculating them now")
        alignment_result = build_pairwise_alignments(sequences)
    pairwise_alignments = alignment_result["pairwise_alignments"]
    min_dist_pair = alignment_result["min_dist_pair"]

    pairwise_alignments = fix_odd_sequence_count(sequences, pairwise_alignments, min_dist_pair)

    # add edges between all sequences weighted by the distances in the distance map
    for i, pairwise_alignments_i in pairwise_alignments.items():
        for j, (seq_i, seq_j, d) in pairwise_alignments_i.items():
            g.add_edge(i, j, weight=(m - d))

    # calculate a max-weight matching
    mate = max_weight_matching(g)

    for i, j in mate:
        if i > j:
            i, j = j, i
        pairwise_alignment = pairwise_alignments[i][j]

        # if f is part of the sequence match, replace it by its sources a and b and
        # include their distance d in the distance result
        if j == len(sequences):
            a, b, d = min_dist_pair
            generalized_sequences.append({
                "group": [pairwise_alignment[0].id, sequences[a].id, sequences[b].id],
                "generalizedSeq": generalize(seq_str(pairwise_alignment[0]), seq_str(pairwise_alignment[1])),
                "distance": pairwise_alignment[2] + d
            })
            generalized_distance = generalized_distance + pairwise_alignment[2] + d
        else:
            generalized_sequences.append({
                "group": [pairwise_alignment[0].id, pairwise_alignment[1].id],
                "generalizedSeq": generalize(seq_str(pairwise_alignment[0]), seq_str(pairwise_alignment[1])),
                "distance": pairwise_alignment[2]
            })
            generalized_distance = generalized_distance + pairwise_alignment[2]

    return generalized_sequences, generalized_distance
