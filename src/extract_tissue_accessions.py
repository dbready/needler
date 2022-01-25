# Extract the identified Liver proteins from the Tissue Atlas paper.
# Use the Uniprot webservice to convert from Ensembl protein ID into uniprot accession

import argparse
import urllib.parse
import urllib.request

import pandas as pd
import tqdm

URL = "https://www.uniprot.org/uploadlists/"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tissuedb", dest="tissuedb")
    parser.add_argument("dst")

    args = parser.parse_args()
    db = pd.read_excel(args.tissuedb, sheet_name="A. Protein copies")

    liver_proteins = db[db["Liver"] > 0]
    # pull the Ensembl Protein ID
    prot_ids = list(set(liver_proteins["Protein ID"]))

    # use the uniprot webservice to convert Ensembl <-> Accession
    # query subset at a time
    BLOCK_SIZE = 50
    rows = []
    for pivot in tqdm.tqdm(range(0, len(prot_ids), BLOCK_SIZE)):
        query_ids = " ".join(prot_ids[pivot : pivot + BLOCK_SIZE])
        params = {
            "from": "ENSEMBL_PRO_ID",
            "to": "ACC",
            "format": "tab",
            "query": query_ids,
        }
        data = urllib.parse.urlencode(params)
        data = data.encode("utf-8")
        request = urllib.request.Request(URL, data)
        with urllib.request.urlopen(request) as handle:
            response = handle.read()
            for i, line in enumerate(response.decode("utf-8").splitlines()):
                if i != 0:
                    ensembl, accession = line.split("\t")
                    rows.append(
                        {
                            "ensembl_protein_id": ensembl.strip(),
                            "accession": accession.strip(),
                        }
                    )

    df = pd.DataFrame(rows)
    df.drop_duplicates(inplace=True)
    df.to_csv(args.dst, index=False, encoding="utf8")


if __name__ == "__main__":
    main()
