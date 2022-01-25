import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pepdb", dest="pepdb")
    parser.add_argument("--filtersrc", dest="filter_src")
    parser.add_argument("dstdb")

    args = parser.parse_args()
    db = pd.read_csv(args.pepdb)
    filter_src = pd.read_csv(args.filter_src)

    accessions = set(filter_src.accession)
    res = db[db.accession.isin(accessions)]
    res.to_csv(args.dstdb, index=False, encoding="utf8", float_format="%.4f")


if __name__ == "__main__":
    main()
