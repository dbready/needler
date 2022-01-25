#!/usr/bin/env python

# Use the z3 constraint solver to find the collection of peptides that maximally utilize a MS gradient.

# Note that we are using a bit of a hack- the python z3 interface does offer an
# Optimizer to which we can apply a timeout function. However, we do an ersatz
# optimizer, we record the best observed model, tell z3 to do better than that
# one, and let it continue until the timeout. Nice property is we can observe 
# how the model is progressing in real-time and save the best model as it is generated.

import argparse
import datetime
import logging
import os
import random
import shutil
import uuid

import pandas as pd
import z3

logger = logging.getLogger()


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dst", dest="dst", required=True, help="csv result filename")
    parser.add_argument("--shuffle", dest="shuffle", action="store_true")
    parser.add_argument("--sortgroups", dest="sort_groups", action="store_true")
    parser.add_argument(
        "--timeout",
        dest="timeout",
        default=-1,
        type=int,
        help="maximum solving time in seconds (-1 unlimited)",
    )
    parser.add_argument(
        "--targets",
        dest="targets_per_cycle",
        type=int,
        default=10,
        help="targets per second (default 10)",
    )
    parser.add_argument(
        "--pepsperprot",
        dest="peps_per_prot",
        type=int,
        default=2,
        help="minimum targets per protein (default 2)",
    )
    parser.add_argument(
        "--rtwidth",
        dest="rt_width",
        type=int,
        help="+/- seconds surrounding the rt column",
    )
    parser.add_argument("--verbose", dest="verbose", action="store_true")
    parser.add_argument("infile")
    return parser.parse_args()


class Needler:
    def __init__(
        self,
        db_src,
        dst,
        rt_width,
        targets_per_cycle,
        timeout=-1,
        peps_per_prot=2,
        shuffle=True,
        sort_groups=False,
    ):
        self.db_src = db_src
        self.dst = dst
        self.rt_width = rt_width
        self.targets_per_cycle = targets_per_cycle
        self.peps_per_prot = peps_per_prot
        self.timeout = timeout
        self.shuffle = shuffle
        self.sort_groups = sort_groups

    def load(self):
        self.df = pd.read_csv(self.db_src)
        initial_rows = self.df.shape[0]
        self.df["protein_id"] = self.df.accession
        self.df["peptide_id"] = self.df.peptide_sequence

        logger.debug(f"Input file unique proteins: {self.df.protein_id.nunique()}")
        logger.debug(f"Input file unique peptides: {self.df.peptide_id.nunique()}")

        self.df["protein_peps"] = self.df.groupby(["protein_id"]).peptide_id.transform(
            lambda x: x.nunique()
        )
        # are there any proteins with fewer than minimum peptides_per_protein? remove them
        bidx = self.df.protein_peps < self.peps_per_prot
        if bidx.sum() > 0:
            acc_count_removed = self.df[bidx].protein_id.nunique()
            self.df = self.df.loc[self.df.protein_peps >= args.peps_per_prot]
            logger.info(
                f"Dropping proteins with insufficient available peptides, proteins removed: {acc_count_removed}"
            )

        self.protein_order = sorted(set(self.df.protein_id))
        self.peptide_order = sorted(set(self.df.peptide_id))
        if self.shuffle:
            # z3 is ~deterministic, some ordering is dependent upon
            # variable name, contents, and/or insert order. Using uuids to help randomize
            accession_guid = self.df.groupby(["protein_id"]).apply(
                lambda x: "{}".format(uuid.uuid4())
            )
            self.df["protein_id"] = self.df.protein_id.replace(accession_guid)
            self.df["peptide_id"] = self.df["peptide_id"].apply(
                lambda x: f"{uuid.uuid4()}"
            )
            logger.debug("Shuffling the protein and peptide ids")
            # re-order dataframe rows
            self.df = self.df.sample(frac=1).reset_index(drop=True)
            # dict insertion ordering preserves and removes duplicates
            self.peptide_order = list(dict.fromkeys(self.df.peptide_id))
            self.protein_order = list(dict.fromkeys(self.df.protein_id))
        if self.sort_groups:
            logger.debug("Shuffling the peptide ids with protein size-ordering")
            # stable sort the proteins, prioritizing those which
            # have fewer pep/protein, maintain sub-ordering
            # ie, all of the proteins with say 2 peptides will be at the top, followed by the 3s, etc
            # does not sub-group the peptides within a protein
            # -invoking this option means titin will allways be at the end
            self.df.sort_values(
                ["protein_peps"], ascending=[True], kind="mergesort", inplace=True
            )
            self.protein_order = list(dict.fromkeys(self.df.protein_id))

        self.df["rt_start"] = (
            self.df.prosit_predicted_rt_seconds - self.rt_width
        ).astype(int)
        self.df["rt_stop"] = (
            self.df.prosit_predicted_rt_seconds + self.rt_width
        ).astype(int)
        bidx = self.df.rt_start < 0
        if bidx.sum() > 0:
            logger.info(
                f"Adjusting negative peptide rt_start values to 0 for peptides: {bidx.sum()}"
            )

        if initial_rows != self.df.shape[0]:
            logger.debug(
                "Filtering occured, remaining unique proteins: {self.df.protein_id.nunique()}"
            )
            logger.debug(
                "Filtering occured, remaining peptides: {self.df.peptide_id.nunique()}"
            )

    def build_model(self):
        self.solver = z3.Solver()
        # enabling these two configurations allows speedups on using pseudo-booleans
        self.solver.set("sat.pb.solver", "solver")
        self.solver.set("sat.cardinality.solver", True)
        if self.timeout > 0:
            logger.debug(f"Setting z3 solver timeout seconds: {self.timeout}")
            # z3 api expects milliseconds
            self.solver.set("timeout", int(self.timeout) * 1000)

        logger.info("Creating z3 Peptide representations")
        self.peptides = {}
        for pep_id in self.peptide_order:
            # implicitly benefiting from those ordered dictionaries
            self.peptides[pep_id] = z3.Bool(pep_id)

        # give z3 a hint that at most XXX peptides should be true in a given solution
        # For example: gradient is 90 seconds long, 3 peptides/cycle, peptide elutes over 30 seconds: can only target 9 peps
        # Alternatively, if using very small windows (eg 1 second) this becomes huge number and want to use
        # the #proteins*targets_per_protein
        # testing showed having an offset improved performance
        PEPTIDE_ALLOWED_OFFSET = 500
        pep_elution_width = 1 + (2 * self.rt_width)
        gradient_length_seconds = self.df.rt_stop.max() - self.df.rt_start.min()
        blocks = (gradient_length_seconds // pep_elution_width) + 1
        max_targetable_blocks = blocks * self.targets_per_cycle
        max_target_from_prots = len(self.protein_order) * self.peps_per_prot
        max_targetable = (
            min(max_targetable_blocks, max_target_from_prots) + PEPTIDE_ALLOWED_OFFSET
        )
        logger.debug(f"Modeling max targetable peptides (including offset) as: {max_targetable}")
        pep_bools = [self.peptides[pid] for pid in self.peptide_order]
        self.solver.add(z3.AtMost(*pep_bools, max_targetable))

        logger.info("Assigning the protein-peptide relationship")
        self.proteins_targeted = []
        for prot_id in self.protein_order:
            protein_peptide_ids = set(self.df[self.df.protein_id == prot_id].peptide_id)
            # PseudoBoolean tricks require returning a tuple of the (Z3.Bool, int) representing value
            pseudo_bools = [(self.peptides[pid], 1) for pid in protein_peptide_ids]
            if args.shuffle:
                random.shuffle(pseudo_bools)
            # collection of peptides associated with this protein must be either
            # - 0 (none selected) or equal to the target value
            self.solver.add(
                z3.Or(
                    z3.PbEq(pseudo_bools, 0), z3.PbEq(pseudo_bools, self.peps_per_prot)
                )
            )
            # this is collection for monitoring if the protein is targeted
            # monitor this collection for overall goal progress
            self.proteins_targeted.append(
                z3.If(z3.PbEq(pseudo_bools, self.peps_per_prot), 1, 0)
            )

        # possible to define larger starting point
        protein_target = 0
        logger.info(f"Setting intial protein target >\t{protein_target}")
        self.solver.add(z3.Sum(self.proteins_targeted) > protein_target)

        # Iterate over each cycle of the gradient (eg seconds)
        # Identify all peptides eluting at this point
        # Sum any targeted peptides and ensure that the number of peptides is
        # less than the maximum allowed targets per cycle
        logger.info("Assigning gradient constraints")
        rt_minimum = max(self.df.rt_start.min(), 0)
        rt_maximum = self.df.rt_stop.max()
        time_steps = list(range(rt_minimum, rt_maximum + 1))
        logger.debug(f"Time steps to model: {len(time_steps)}")
        if self.shuffle:
            random.shuffle(time_steps)
        for i in time_steps:
            # note this is inclusive on both ends of the time range
            bidx = (self.df.rt_start <= i) & (self.df.rt_stop >= i)
            if bidx.sum() > 0:
                pep_ids = set(self.df[bidx].peptide_id)
                time_slot_peptides = [self.peptides[pid] for pid in pep_ids]
                if args.shuffle:
                    random.shuffle(time_slot_peptides)
                self.solver.add(z3.AtMost(*time_slot_peptides, self.targets_per_cycle))

    def save(self):
        # write to temp name in case there is crash
        temp_name = self.dst + ".temp"
        parent_dir = os.path.dirname(temp_name)
        if not os.path.exists(parent_dir):
            os.makedirs(parent_dir, exist_ok=True)
        # only save the valid peptides & proteins
        bidx = (self.df["peptide_in_model"]) & (self.df["protein_is_targeted"])
        columns = [
            "accession",
            "isoform",
            "name",
            "gene",
            "peptide_sequence",
            "peptide",
            "prosit_predicted_rt_seconds",
            "rt_start",
            "rt_stop",
        ]
        self.df[bidx][columns].to_csv(temp_name, index=False, encoding="utf8")
        shutil.move(temp_name, self.dst)

    def optimize(self):
        logger.info("Begin model evaluation")
        model = None
        # z3 timeout functionality does not work with this re-evaluation mechanism, have to
        # re-assign the timeout between solutions
        start_dt = datetime.datetime.now()
        target_dt = start_dt + datetime.timedelta(seconds=self.timeout)
        best_ever_result = 0
        while True:
            # sat, unsat, unknown
            result = self.solver.check()

            # reset the timeout budget
            if self.timeout > 0:
                remainder_seconds = (
                    target_dt - datetime.datetime.now()
                ).total_seconds()
                if remainder_seconds < 0:
                    # using 1ms to trigger it if we somehow managed to slip over the 0 threshold
                    # because we want to do the rest of our cleanup actions first
                    self.solver.set("timeout", 1)
                else:
                    self.solver.set("timeout", int(remainder_seconds) * 1000)

            if result == z3.sat:
                model = self.solver.model()
                used_pep_ids = set()
                for peptide_id, pep_bool in self.peptides.items():
                    was_used = z3.is_true(model[pep_bool])
                    if was_used:
                        used_pep_ids.add(peptide_id)
                self.df["peptide_in_model"] = self.df.peptide_id.isin(used_pep_ids)
                # peptide is used...does it define a fully valid protein with min peptides_per_prot?
                acc_counts = (
                    self.df[self.df.peptide_in_model]
                    .groupby(["protein_id"])
                    .peptide_id.nunique()
                )
                protein_targets = acc_counts[acc_counts == self.peps_per_prot].index
                self.df["protein_is_targeted"] = self.df.protein_id.isin(
                    protein_targets
                )
                protein_target = len(set(protein_targets))
                logger.debug(
                    "Created model targeting proteins: {}\tpeptides: {}".format(
                        protein_target, len(used_pep_ids)
                    )
                )
                if protein_target > best_ever_result:
                    best_ever_result = protein_target
                    # save best model to date
                    self.save()
                    ## set harder target for the solver
                    self.solver.add(z3.Sum(self.proteins_targeted) > protein_target)
                else:
                    logger.warning(
                        f"Z3 loop detected, best seen was:\t{best_ever_result}"
                    )
                    break
                if protein_target == self.df.protein_id.nunique():
                    logger.info("SOLVED! Model targets all provided proteins")
                    break
            elif result == z3.unsat:
                logger.info(
                    "Impossible to generate a satisfactory model. Relax the constraints."
                )
                break
            elif result == z3.unknown:  # 'unknown' when z3 timeout hit or CTRL-C
                if remainder_seconds < 0:
                    logger.error("Model timeout reached")
                else:
                    logger.error(f"Model generation was halted, status '{result}")
                break


def main(*args, **kwargs):

    needle = Needler(*args, **kwargs)

    # load datafrane and do munging
    needle.load()

    # build z3 model representation
    needle.build_model()

    # begin solving
    needle.optimize()


if __name__ == "__main__":
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s\t%(levelname)-8s\t%(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    args = get_args()
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    logger.debug(args)

    main(
        db_src=args.infile,
        dst=args.dst,
        rt_width=args.rt_width,
        targets_per_cycle=args.targets_per_cycle,
        timeout=args.timeout,
        peps_per_prot=args.peps_per_prot,
        shuffle=args.shuffle,
        sort_groups=args.sort_groups,
    )
