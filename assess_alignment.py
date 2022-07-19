#!/usr/bin/env python

import argparse
import json
import sys
from dataclasses import dataclass
from typing import TextIO


@dataclass
class Sample:
    name: str
    reference: str
    mismatches: int
    quality: str


def load_json(json_file: TextIO) -> dict:
    data = json.load(json_file)
    if "msa" not in data:
        raise ValueError("MSA missing from JSON, cannot proceed")
    if "gappedTraces" not in data:
        raise ValueError("gappedTraces missing from JSON, cannot proceed")
    return data


def analyse_mismatches(
    json_file: TextIO,
    offset: int,
    length: int,
    vp1only: bool = True,
    sec_is_conflict: bool = False,
) -> list:
    data = load_json(json_file)
    msas = [al for al in data["msa"] if not al["reference"]]
    reference = None
    for al in data["msa"]:
        if al["reference"]:
            if reference is None:
                reference = al
            else:
                sys.exit(
                    "more than one reference found in JSON MSA list, cannot proceed"
                )
    min_start = min([int(al["leadingGaps"]) for al in msas])
    max_end = max([int(al["leadingGaps"]) + len(al["align"]) for al in msas])
    base_state = ["n"] * len(reference["align"])
    mismatch_bases = {}
    for i, base in enumerate(reference["align"]):
        for k, al in enumerate(msas):
            leading_gaps = int(al["leadingGaps"])
            align_len = len(al["align"])
            if leading_gaps < i and (leading_gaps + align_len) > i:
                vp1pos = i - offset
                if vp1only and vp1pos < 0 or vp1pos > length:
                    # skip positions outside of vp1 gene region
                    continue
                al_base = al["align"][i - leading_gaps]
                has_secondary_basecall = False
                if sec_is_conflict:
                    gappedTrace = data["gappedTraces"][k]
                    pos = i - int(gappedTrace["leadingGaps"])
                    # print(len(gappedTrace['basecallPos']), pos, k, len(gappedTrace['basecalls']), gappedTrace['basecallPos'][pos])
                    basecall_str = gappedTrace["basecalls"][
                        str(gappedTrace["basecallPos"][pos])
                    ]
                    if "|" in basecall_str:
                        has_secondary_basecall = True
                        # set this position to conflicted
                        base_state[i] = "C"
                if al_base != base:
                    # let's deal with all the cases where the base state doesn't match the reference
                    if base_state[i] == "G":
                        # the base state was G (a trace matches reference) and now we see a mismatch
                        base_state[i] = "C"
                    elif base_state[i] == "C":
                        # already marked as conflicting - a mismatch doesn't change that
                        pass
                    elif base_state[i] == "n" or base_state[i] == "M":
                        # we never saw this before or its already marked as a mismatch
                        base_state[i] = "M"
                        mismatch_bases[i] = al_base
                    else:
                        sys.exit("unexpected base state: " + base_state[i])
                else:
                    if base_state[i] == "G" or base_state[i] == "n":
                        # we saw this before and got a match or
                        # we never saw this before
                        base_state[i] = "G"
                    elif base_state[i] == "M":
                        # we saw this before but it was a mismatch - mark this as a conflict
                        base_state[i] = "C"
                        if i in mismatch_bases:
                            del mismatch_bases[i]
                    elif base_state[i] == "C":
                        # we have seen a conflict here before
                        pass
                    else:
                        sys.exit("unexpected base_state: " + base_state[i])
    conflicts = base_state.count("C")
    matches = base_state.count("G")
    mismatches = base_state.count("M")
    mismatch_list = []
    for i, state in enumerate(base_state):
        # i is in zero-based genome coordinates
        if state == "M":
            # for mismatch store [pos_in_genome, pos_in_vp1, reference_base, sequenced_base]
            mismatch_list.append(
                [i, i - offset, reference["align"][i], mismatch_bases[i]]
            )
    return [conflicts, matches, mismatches, mismatch_list]


def analyse_trace_quality(json_file: TextIO) -> float:
    data = load_json(json_file)

    traces = data["gappedTraces"]
    overall_avg = 0
    for trace in traces:
        start = min(trace["basecallPos"])
        end = max(trace["basecallPos"])
        call_quality = {}
        avg_ratio = 0
        for base in ("A", "C", "G", "T"):
            calls = trace["peak" + base][start : end + 1]
            min_call = min(calls)
            max_call = max(calls)
            avg_call = sum(calls) / len(calls)
            ratio = max_call / avg_call
            call_quality["avg" + base] = avg_call
            call_quality["min" + base] = min_call
            call_quality["max" + base] = max_call
            call_quality["ratio" + base] = ratio
            avg_ratio += ratio
        avg_ratio = avg_ratio / 4
        overall_avg += avg_ratio
    overall_avg = overall_avg / len(traces)
    return overall_avg


def comma_split(args: str) -> list[str]:
    return args.split(",")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_filename", help="Path to output file")
    parser.add_argument("--sample_name", help="Name of sample being analysed")
    parser.add_argument(
        "--dataset_names", type=comma_split, help="Comma separated names for datasets"
    )
    parser.add_argument("--datasets", nargs="+")
    args = parser.parse_args()

    offsets = {
        "poliovirus1sabin": 2480,
        "poliovirus2sabin": 2482,
        "poliovirus3sabin": 2477,
    }

    lengths = {
        "poliovirus1sabin": 906,
        "poliovirus2sabin": 903,
        "poliovirus3sabin": 900,
    }

    min_mismatches = None
    for file_index, json_filename in enumerate(args.datasets):
        dataset_name = args.dataset_names[file_index].replace(
            ".json", ""
        )  # take the name but remove any json suffix
        offset = offsets[dataset_name]
        length = lengths[dataset_name]
        (conflicts, matches, mismatches, mismatch_list) = analyse_mismatches(
            open(json_filename), offset, length
        )
        # analyse_mismatches(json_filename, True)
        quality = analyse_trace_quality(open(json_filename))
        if min_mismatches is None or mismatches < min_mismatches:
            min_mismatches = mismatches
            best_match_mismatch_list = mismatch_list
            best_match_quality = quality

    info = {
        "sample_name": args.sample_name,
        "best_reference": dataset_name,
        "mismatches": min_mismatches,
        "mismatch_list": best_match_mismatch_list,
        "quality": best_match_quality,
    }
    json.dump(info, open(args.output_filename, "w"))
