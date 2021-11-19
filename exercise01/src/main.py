import argparse
import json
import pickle

from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline

from .ga import mwm_ga
from .mwm import mwm_exact, build_pairwise_alignments
from .dnala import dnala

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generalize DNA sequences with DNALA, exact MWM, or approximate MWM')
    parser.add_argument('INFILE', help='the FASTA file that holds the DNA sequences to be generalized')
    parser.add_argument('-a', '--algorithm', default='DNALA', help='the algorithm to apply; DNALA | MWM | GA | PWA_ONLY')
    parser.add_argument('-al', '--alignments', nargs='?', help='a file containing pre-computed alignments for the given INFILE')
    parser.add_argument('-o', '--outfile', nargs='?', help='the file to write the results to')

    args = parser.parse_args()
    # sequences = list(SeqIO.parse(args.INFILE, "fasta"))

    if args.algorithm == "DNALA":
        if args.alignments is None:
            print(f"---- Aligning sequences in file {args.INFILE} with CLUSTALW ----")
            clustal_cline = ClustalwCommandline("clustalw2", infile=args.INFILE, outorder="INPUT")
            clustal_cline()
            alignment_filename = args.INFILE.rpartition(".")[0] + ".aln"
            aligned_sequences = list(AlignIO.read(alignment_filename, "clustal"))
        else:
            aligned_sequences = list(AlignIO.read(args.alignments, "clustal"))

        print(f"---- Executing DNALA algorithm on aligned sequences ----")
        gen_sequences, total_distance = dnala(aligned_sequences)

        print("---- Results ----")
        print("Total generalized distance: ")
        print(total_distance)
        print("Generalized Sequences: ")
        print(json.dumps(gen_sequences, indent=2))

        if args.outfile is not None:
            with open(args.outfile, "w") as f:
                f.write(json.dumps(gen_sequences, indent=2))

    elif args.algorithm == "MWM":
        sequences = list(SeqIO.parse(args.INFILE, "fasta"))

        alignment_result = None
        if args.alignments is not None:
            with open(args.alignments, "rb") as f:
                alignment_result = pickle.load(f)

        print(f"---- Executing MWM algorithm on aligned sequences ----")
        gen_sequences, total_distance = mwm_exact(sequences, alignment_result)

        print("---- Results ----")
        print("Total generalized distance: ")
        print(total_distance)
        print("Generalized Sequences: ")
        print(json.dumps(gen_sequences, indent=2))

        if args.outfile is not None:
            with open(args.outfile, "w") as f:
                f.write(json.dumps(gen_sequences, indent=2))

    elif args.algorithm == "PWA_ONLY":
        sequences = list(SeqIO.parse(args.INFILE, "fasta"))
        alignment_result = build_pairwise_alignments(sequences)

        with open(args.outfile, "wb") as f:
            pickle.dump(alignment_result, f)

    elif args.algorithm == "GA":
        sequences = list(SeqIO.parse(args.INFILE, "fasta"))

        alignment_result = None
        if args.alignments is not None:
            with open(args.alignments, "rb") as f:
                alignment_result = pickle.load(f)

        print(f"---- Executing GA algorithm on aligned sequences ----")
        gen_sequences, total_distance = mwm_ga(sequences, alignment_result)

        print("---- Results ----")
        print("Total generalized distance: ")
        print(total_distance)
        print("Generalized Sequences: ")
        print(json.dumps(gen_sequences, indent=2))

        if args.outfile is not None:
            with open(args.outfile, "w") as f:
                f.write(json.dumps(gen_sequences, indent=2))
