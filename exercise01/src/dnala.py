import math

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


def distance_matrix(sequences: [str]) -> [[int]]:
    return [[dist(seq_a, seq_b) for seq_b in sequences] for seq_a in sequences]


def dnala(sequences: [str]) -> [str]:
    dist_matrix = distance_matrix(sequences)
    generalized_sequences = []
    remaining_sequences = list(range(len(sequences)))

    def closest_sequences(i):
        min_dist = min(dist_matrix[i][j] if i != j else math.inf for j in remaining_sequences)
        return [j for j in remaining_sequences if dist_matrix[i][j] == min_dist and j != i]

    while len(remaining_sequences) > 0:
        print(remaining_sequences)
        if(len(remaining_sequences)) == 3:
            [x, y, z] = remaining_sequences
            generalized_sequences.append(generalize(generalize(sequences[x], sequences[y]), sequences[z]))
            remaining_sequences = []

        for s in remaining_sequences:
            candidates = closest_sequences(s)
            match_found = False
            for c in candidates:
                reverse_candidates = closest_sequences(c)
                if s in reverse_candidates:
                    print(f"Matched {sequences[s]} to {sequences[c]}: {generalize(sequences[s], sequences[c])}")
                    generalized_sequences.append(generalize(sequences[s], sequences[c]))
                    remaining_sequences.remove(s)
                    remaining_sequences.remove(c)
                    match_found = True
                    break
            if match_found:
                break

    return generalized_sequences
