import argparse
import json

from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline

from .mwm import mwm_exact
from .dnala import dnala

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generalize DNA sequences with DNALA, exact MWM, or approximate MWM')
    parser.add_argument('INFILE', help='the FASTA file that holds the DNA sequences to be generalized')
    parser.add_argument('-a', '--algorithm', default='DNALA', help='the algorithm to apply; DNALA | MWM_EXACT | MWM_GA')

    args = parser.parse_args()
    # sequences = list(SeqIO.parse(args.INFILE, "fasta"))

    if args.algorithm == "DNALA":
        print(f"---- Aligning sequences in file {args.INFILE} with CLUSTALW ----")
        clustal_cline = ClustalwCommandline("clustalw2", infile=args.INFILE, outorder="INPUT")
        clustal_cline()
        alignment_filename = args.INFILE.rpartition(".")[0] + ".aln"
        aligned_sequences = list(AlignIO.read(alignment_filename, "clustal"))

        print(f"---- Executing DNALA algorithm on aligned sequences ----")
        gen_sequences, total_distance = dnala(aligned_sequences)

        print("---- Results ----")
        print("Total generalized distance: ")
        print(total_distance)
        print("Generalized Sequences: ")
        print(json.dumps(gen_sequences, indent=2))

    elif args.algorithm == "MWM_EXACT":
        sequences = list(SeqIO.parse(args.INFILE, "fasta"))

        print(f"---- Executing MWM EXACT algorithm on aligned sequences ----")
        gen_sequences, total_distance = mwm_exact(sequences)

        print("---- Results ----")
        print("Total generalized distance: ")
        print(total_distance)
        print("Generalized Sequences: ")
        print(json.dumps(gen_sequences, indent=2))
