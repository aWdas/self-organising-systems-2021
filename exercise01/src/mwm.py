from Bio import Align, SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
import networkx as nx
from networkx.algorithms.matching import max_weight_matching

from .util import dist, seq_str, generalize


def pairwise_align(seq_a: SeqRecord, seq_b: SeqRecord):
    # this is slow af, but I haven't found a way to make the pairwise2 matcher behave the same as CLUSTAL
    # (it always introduces gaps for sequences of the same length)
    print(f"Pairwise aligning sequence {seq_a.id} to {seq_b.id}")
    SeqIO.write([seq_a, seq_b], "resources/scratch.fasta", "fasta")
    clustal_cline = ClustalwCommandline("clustalw2", infile="resources/scratch.fasta", outorder="INPUT")
    clustal_cline()
    aligned_sequences = list(AlignIO.read("resources/scratch.aln", "clustal"))
    return aligned_sequences[0], aligned_sequences[1]


def mwm_exact(sequences: [SeqRecord]):
    g = nx.Graph()
    m = 10000000000
    generalized_sequences = []
    generalized_distance = 0

    # build a map of the pairwise distances between the sequences,
    # and keep track of the pair with the minimum distance
    pairwise_alignments = {}
    min_dist_pair = None
    for i in range(len(sequences)):
        pairwise_alignments_i = {}
        for j in range(i+1, len(sequences)):
            aligned_seq_i, aligned_seq_j = pairwise_align(sequences[i], sequences[j])
            d = dist(seq_str(aligned_seq_i), seq_str(aligned_seq_j))
            pairwise_alignments_i[j] = (aligned_seq_i, aligned_seq_j, d)
            if min_dist_pair is None or min_dist_pair[2] > d:
                min_dist_pair = (i, j, d)
        pairwise_alignments[i] = pairwise_alignments_i

    # if there's an odd number of sequences, replace the min distance pair by its generalization f in the distance map
    if len(sequences) % 2 == 1:
        a, b, d = min_dist_pair
        print(f"Min-Distance Pair: {min_dist_pair}")
        # remove sequences a and b from the distance map
        pairwise_alignments.pop(a, None)
        pairwise_alignments.pop(b, None)
        for pairwise_alignments_i in pairwise_alignments.values():
            pairwise_alignments_i.pop(a, None)
            pairwise_alignments_i.pop(b, None)

        # add f with index "len(sequences)" to the distance map, calculate the pairwise distances to the other seqs
        seq_f = SeqRecord(Seq(generalize(seq_str(sequences[a]), seq_str(sequences[b]))), id="f")
        for (i, pairwise_alignments_i) in pairwise_alignments.items():
            aligned_seq_i, aligned_seq_f = pairwise_align(sequences[i], seq_f)
            d = dist(seq_str(aligned_seq_i), seq_str(aligned_seq_f))
            pairwise_alignments_i[len(sequences)] = (aligned_seq_i,  aligned_seq_f, d)

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
