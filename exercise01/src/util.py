from Bio.SeqRecord import SeqRecord

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
