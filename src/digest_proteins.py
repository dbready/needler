import argparse
from collections import defaultdict
import os

import pandas as pd
from pyteomics import parser as digest
from pyteomics import mass


def composition_to_string(x):
    pieces = []
    for atom in sorted(x.keys()):
        pieces.append("{}{}".format(atom, x[atom]))
    return "".join(pieces)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--proteins", dest="proteins")
    parser.add_argument("--min_length", dest="min_length", type=int, default=5)
    parser.add_argument("--max_length", dest="max_length", type=int, default=30)
    parser.add_argument("--min_peptides", dest="min_peptides", type=int, default=2)
    parser.add_argument("--unfiltered_filename", dest="unfiltered_filename")
    parser.add_argument(
        "peptides",
    )

    args = parser.parse_args()
    return args


def main():
    args = get_args()

    proteins = pd.read_csv(args.proteins)
    dst_filename = os.path.normpath(args.peptides)

    # create dictionary keyed on peptide sequence with list of protein accessions
    peptides = defaultdict(list)
    for _, row in proteins.iterrows():
        accession = row["accession"]
        sequence = row["sequence"]
        peps = digest.cleave(
            sequence,
            digest.expasy_rules["trypsin"],
            missed_cleavages=0,
            min_length=args.min_length,
        )
        for pep in peps:
            # some proteins have 'X' and 'U' (selenocystine) amino acids listed
            # have no mechanism to handle those currently
            if len(pep) <= args.max_length and "X" not in pep and "U" not in pep:
                peptides[pep].append(accession)
    # identify the peptides represented by single protein accession
    distinct = []
    for peptide, accessions in peptides.items():
        if len(accessions) == 1:
            distinct.append({"accession": accessions[0], "peptide_sequence": peptide})
    distinct_peps = pd.DataFrame(distinct)
    proteins.drop(["sequence", "taxonomy_id"], axis=1, inplace=True)
    df = pd.merge(proteins, distinct_peps, how="inner", on="accession")
    # all peptides are carbamidomethylated
    df["peptide"] = df.peptide_sequence.str.replace("C", "cC")
    # build the modded cys definition
    aa_comp = dict(mass.std_aa_comp)
    aa_comp["cC"] = aa_comp["C"] + mass.Composition("C2H3NO")
    compositions = df.peptide.apply(
        lambda x: mass.Composition(sequence=x, aa_comp=aa_comp)
    )
    df["formula"] = compositions.apply(composition_to_string)
    df["mw"] = compositions.apply(
        lambda x: mass.calculate_mass(x, average=False, charge=0)
    )

    df.sort_values(["gene", "accession", "isoform", "peptide_sequence"], inplace=True)
    if args.unfiltered_filename is not None:
        df.to_csv(args.unfiltered_filename, index=False, encoding="utf8")

    # drop Met containing peptides, unreliable targets not worth considering
    df = df.loc[~df.peptide_sequence.str.contains("M")]
    # find proteins represented by >= threshold peptides
    acc_counts = df.groupby(["accession"]).peptide_sequence.nunique()
    accessions_gte_peps = list(acc_counts[acc_counts >= args.min_peptides].index)

    df[df.accession.isin(accessions_gte_peps)].to_csv(
        dst_filename, index=False, encoding="utf8", float_format="%.4f"
    )


if __name__ == "__main__":
    main()
