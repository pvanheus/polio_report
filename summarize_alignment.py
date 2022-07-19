#!/usr/bin/env python

import argparse
import json

def comma_split(args: str) -> list[str]:
    return args.split(",")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary_output_filename", help="Path to summary output file")
    parser.add_argument("--variant_list_outputs", nargs="+")
    parser.add_argument("--datasets", nargs="+")
    args = parser.parse_args()

    summary_output = open(args.summary_output_filename, "w")
    for i, json_filename in enumerate(args.datasets):
        alignment_summary = json.load(open(json_filename))
        mismatch_list = ";".join(
            [
                ":".join([str(part) for part in el[1:]])
                for el in alignment_summary["mismatch_list"]
            ]
        )
        print(
            alignment_summary["sample_name"],
            alignment_summary["best_reference"],
            alignment_summary["mismatches"],
            round(alignment_summary["quality"], 2),
            mismatch_list,
            sep="\t",
            file=summary_output,
        )
        variant_list_output = open(args.variant_list_outputs[i], "w")
        print('genome pos', 'VP1 pos', 'ref', 'sequence', sep='\t', file=variant_list_output)
        for variant in alignment_summary["mismatch_list"]:
            print("\t".join([str(el) for el in variant]), file=variant_list_output)
        variant_list_output.close()
