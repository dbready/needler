# Extract .csv representation of proteins from .fasta

import argparse
import os

import numpy as np
import pandas as pd
from pyteomics import fasta

SPECIES = {
    "OX=9606",  # human
    # "OX=10090",  # mouse
    # "OX=10116",  # rat
    # "OX=83333",  # ecoli
    # "OX=559292",  # yeast
}


def gen_entries(filename):
    with fasta.read(filename) as handle:
        for header, seq in handle:
            for species_name in SPECIES:
                if species_name in header:
                    yield (header, seq)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", dest="fasta")
    parser.add_argument(
        "dst",
    )

    args = parser.parse_args()

    fasta_filename = os.path.normpath(os.path.abspath(args.fasta))
    dst_filename = os.path.normpath(args.dst)

    proteins = gen_entries(fasta_filename)
    df = pd.DataFrame(list(proteins), columns=("header", "sequence"))
    tokens = df.header.str.split("|", 2)
    df["accession"] = tokens.str[1].str.replace("-1$", "", regex=False)
    # if isoform is not specified, is implicit '1' - make it explicit
    df["isoform"] = tokens.str[1].str.split("-").str[1].replace(np.nan, 1).astype(int)
    df["name"] = tokens.str[2].str.split(" ", 1).str[0]
    df["description"] = tokens.str[2].str.split(" ", 1).str[-1]
    df["taxonomy_id"] = (
        df["description"].str.extract("OX=([\S]+)")[0].str.strip().astype(int)
    )
    df["gene"] = df["description"].str.extract("GN=([\S]+)")[0].str.strip()
    df.drop(["header"], axis=1, inplace=True)
    df.sort_values(["taxonomy_id", "gene", "name", "isoform"], inplace=True)
    col_order = ["accession", "isoform", "name", "gene", "taxonomy_id", "sequence"]
    df[col_order].to_csv(dst_filename, index=False, encoding="utf8")


if __name__ == "__main__":
    main()
