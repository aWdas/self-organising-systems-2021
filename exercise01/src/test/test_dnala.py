from ..dnala import generalize, dist, dnala

def test_generalize():
    assert generalize("A", "A") == "A"
    assert generalize("N", "N") == "N"
    assert generalize("A", "G") == "R"
    assert generalize("A", "C") == "M"
    assert generalize("A", "R") == "R"
    assert generalize("C", "R") == "V"
    assert generalize("A", "-") == "N"
    assert generalize("A", "B") == "N"
    assert generalize("A", "H") == "H"

def test_dist():
    assert dist("AAA", "AAA") == 0
    assert dist("ACC", "CAA") == 6
    assert dist("A", "N") == 3
    assert dist("T", "B") == 2
    assert dist("A", "M") == 1
    assert dist("A", "Y") == 3

def test_dnala():
    test_sequences = [
        "DATB",
        "AATC",
        "TACC",
        "AACC",
        "TAAC",
        "GTGC"
    ]
    assert dnala(test_sequences) == ["AAYC", "TAMC", "DWKB"]

def test_dnala_uneven():
    test_sequences = [
        "DATB",
        "AATC",
        "TACC",
        "AACC",
        "----",
        "TAAC",
        "GTGC"
    ]
    assert dnala(test_sequences) == ["AAYC", "TAMC", "NNNN"]