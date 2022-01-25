import argparse

import numpy as np
import pandas as pd
from scipy import stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--peptides", dest="peptides")
    parser.add_argument("--procal", dest="procal")
    parser.add_argument("--prosit", dest="prosit")
    parser.add_argument("dstdb")

    args = parser.parse_args()
    peptides = pd.read_csv(args.peptides)
    procal = pd.read_csv(args.procal)
    prosit = pd.read_csv(args.prosit)  # prost predicted irts

    # add irts to the emperical procal rts
    procal = pd.merge(procal, prosit, how="left", on="peptide_sequence")

    # calculate the linear fit between the Prosit predicted iRT values
    # and the emperically observed PROCAL peptide RTs
    rt = procal.retention_time_minutes
    irt = procal.prosit_predicted_irt
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(irt, rt)
    print(
        f"slope: {slope}\tintercept: {intercept}\trvalue: {rvalue}\tpvalue: {pvalue}\tstderr:\t{stderr}"
    )

    # align the predicted iRT values for supplied peptides
    res = pd.merge(peptides, prosit, how="left", on="peptide_sequence")
    # convert the iRT<->RT
    rt_minutes = res.prosit_predicted_irt * slope + intercept
    res["prosit_predicted_rt_seconds"] = np.round(rt_minutes * 60.0, 0).astype(int)
    res.to_csv(args.dstdb, index=False, encoding="utf8", float_format="%.4f")


if __name__ == "__main__":
    main()
