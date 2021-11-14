from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import format_alignment

adjacency_list = {
    "A": {"R", "W", "M"},
    "G": {"R", "K", "S"},
    "T": {"W", "K", "Y"},
    "C": {"M", "S", "Y"},
    "R": {"D", "V"},
    "W": {"D", "H"},
    "M": {"V",  "H"},
    "K": {"D", "B"},
    "S": {"V", "B"},
    "Y": {"H", "B"},
    "D": {"N"},
    "V": {"N"},
    "H": {"N"},
    "B": {"N"},
    "-": {"N"},
    "N": {"N"}
}


def g(symbol_a: str, symbol_b: str) -> str:
    def g_internal(symbols_a, symbols_b):
        if len(symbols_a.intersection(symbols_b)) > 0:
            return next(iter(symbols_a.intersection(symbols_b)))
        else:
            with_parents_a = symbols_a.union(*[adjacency_list[k] for k in symbols_a])
            with_parents_b = symbols_b.union(*[adjacency_list[k] for k in symbols_b])

            if len(with_parents_a.intersection(symbols_b)) > 0:
                return next(iter(with_parents_a.intersection(symbols_b)))
            elif len(with_parents_b.intersection(symbols_a)) > 0:
                return next(iter(with_parents_b.intersection(symbols_a)))
            else:
                return g_internal(with_parents_a, with_parents_b)

    return g_internal({symbol_a},  {symbol_b})


def generalize(seq_a: str, seq_b: str) -> str:
    return "".join(g(symbol_a, symbol_b) for symbol_a, symbol_b in zip(seq_a, seq_b))


def level(symbol: str) -> int:
    if symbol in ["A", "G", "T", "C"]:
        return 1
    elif symbol in ["R", "W", "M", "K", "S", "Y"]:
        return 2
    elif symbol in ["D", "V", "H", "B", "-"]:
        return 3
    else:
        return 4


def d(symbol_a: str, symbol_b: str) -> int:
    return 2*level(g(symbol_a, symbol_b)) - level(symbol_a) - level(symbol_b)


def dist(seq_a: str, seq_b: str) -> int:
    return sum(d(symbol_a, symbol_b) for symbol_a, symbol_b in zip(seq_a, seq_b))


def seq_str(sequence: SeqRecord) -> str:
    return str(sequence.seq)

def pairwise_align(seq_a: SeqRecord, seq_b: SeqRecord):
    # this is slow af, but I haven't found a way to make the pairwise2 matcher behave the same as CLUSTAL
    # (it always introduces gaps for sequences of the same length)
    print(f"Pairwise aligning sequence {seq_a.id} to {seq_b.id}")
    # penalize gaps like the default parameters of CLUSTALW according to https://www.genome.jp/tools-bin/clustalw
    alignments = pairwise2.align.globalxs(seq_str(seq_a), seq_str(seq_b), -10.0, -0.1)
    print(format_alignment(*alignments[0]))

    return SeqRecord(Seq(alignments[0].seqA), id=seq_a.id), SeqRecord(Seq(alignments[0].seqB), id=seq_b.id)


def build_pairwise_alignments(sequences: [SeqRecord]):
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

    return {
        "pairwise_alignments": pairwise_alignments,
        "min_dist_pair": min_dist_pair
    }

def fix_odd_sequence_count(sequences, pairwise_alignments, min_dist_pair):
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
            fd = dist(seq_str(aligned_seq_i), seq_str(aligned_seq_f))
            pairwise_alignments_i[len(sequences)] = (aligned_seq_i,  aligned_seq_f, fd)

    return pairwise_alignments

