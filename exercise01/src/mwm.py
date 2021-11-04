from Bio import Align, SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
import networkx as nx
from networkx.algorithms.matching import max_weight_matching

from .util import dist, seq_str, generalize


def pairwise_align(seq_a: SeqRecord, seq_b: SeqRecord):
    # this is slow af, but I haven't found a way to make the pairwise2 matcher behave the same as CLUSTAL
    # (it always introduces gaps for sequences of the same length)
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

    pairwise_alignments = {}

    for i in range(len(sequences)):
        pairwise_alignments_i = {}
        for j in range(i+1, len(sequences)):
            aligned_seq_i, aligned_seq_j = pairwise_align(sequences[i], sequences[j])
            d = dist(seq_str(aligned_seq_i), seq_str(aligned_seq_j))
            pairwise_alignments_i[j] = (aligned_seq_i, aligned_seq_j, d)
            g.add_edge(i, j, weight=(m - d))
        pairwise_alignments[i] = pairwise_alignments_i

    mate = max_weight_matching(g)
    print(mate)
    for i, j in mate:
        if i > j:
            i, j = j, i
        pairwise_alignment = pairwise_alignments[i][j]
        generalized_sequences.append(generalize(seq_str(pairwise_alignment[0]), seq_str(pairwise_alignment[1])))
        generalized_distance = generalized_distance + pairwise_alignment[2]

    return generalized_sequences, generalized_distance
