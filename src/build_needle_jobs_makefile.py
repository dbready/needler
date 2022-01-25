#!/usr/bin/env python

import argparse
from collections import defaultdict
import os

TISSUES = {
    "proteome": "build/pepdb_rts_60min_gradient.csv",
    "liver": "build/liver_rts_60min_gradient.csv",
    "kinase": "build/kinase_rts_60min_gradient.csv",
}


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rts", dest="rts", required=True)
    parser.add_argument("--targets", dest="targets", required=True)
    parser.add_argument(
        "--replicates",
        dest="replicates",
        default=1,
        type=int,
    )
    parser.add_argument("--solve_seconds", default=86_400, dest="solve_seconds", type=int)
    parser.add_argument("--dst_dir", default="build", dest="dst_dir")
    parser.add_argument("dst")
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    rts = args.rts.split(",")
    targets = args.targets.split(",")

    make_targets = defaultdict(list)
    for tissue_name, tissue_db in TISSUES.items():
        for rep_num in range(1, args.replicates + 1):
            for target in targets:
                for rt in rts:
                    filename = os.path.join(
                        args.dst_dir,
                        tissue_name,
                        f"{tissue_name}_TARG{target.zfill(5)}_RT{rt.zfill(3)}_R{rep_num}.csv",
                    )
                    # even without the --verbose flag, might be a bit chatty for running in parallel
                    rule = f"${{PYTHON}} src/needler.py --shuffle --sortgroups --verbose --timeout={args.solve_seconds} --targets={target} --rtwidth={rt} --dst={filename} {tissue_db}"
                    make_targets[tissue_name].append((filename, rule))

    with open(args.dst, "w") as handle:
        # write the top level dependency list for "ease" of human review
        # unreadable, but technically there
        for tissue_name, sub_targets in make_targets.items():
            sub_target_str = " ".join((targ for targ, _ in sub_targets))
            handle.write(f"{tissue_name} : {sub_target_str}\n\n")
        # now write out the actual rules
        for dependencies in make_targets.values():
            for target, rule in dependencies:
                handle.write(f"{target} :\n")
                handle.write(f"\t{rule}\n\n")


if __name__ == "__main__":
    main()
