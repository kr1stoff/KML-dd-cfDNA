import sys
from pathlib import Path
import json
import pandas as pd
import numpy as np

sys.stderr = open(snakemake.log[0], "w")


def fastp_all_samples_qc(files_fastp_json, out_tsv):
    title = ["Sample", "RawReads", "RawBases", "CleanReads", "CleanBases", "RawQ20",
             "RawQ30", "CleanQ20", "CleanQ30", "CleanAverageLength", "GC", "DuplicationRate"]
    df = pd.DataFrame(columns=title)
    for js_path in files_fastp_json:
        js_data = json.loads(open(js_path, "r").read())
        sample = Path(js_path).stem
        mean_lengths = np.array(
            [v for k, v in js_data["summary"]["after_filtering"].items() if k.endswith("mean_length")])
        out = [
            sample,
            js_data["summary"]["before_filtering"]["total_reads"],
            js_data["summary"]["before_filtering"]["total_bases"],
            js_data["summary"]["after_filtering"]["total_reads"],
            js_data["summary"]["after_filtering"]["total_bases"],
            js_data["summary"]["before_filtering"]["q20_rate"],
            js_data["summary"]["before_filtering"]["q30_rate"],
            js_data["summary"]["after_filtering"]["q20_rate"],
            js_data["summary"]["after_filtering"]["q30_rate"],
            int(mean_lengths.mean()),
            js_data["summary"]["after_filtering"]["gc_content"],
            js_data["duplication"]["rate"],
        ]
        df.loc[len(df)] = out
    df.to_csv(out_tsv, index=False, sep="\t")


if __name__ == "__main__":
    fastp_all_samples_qc(snakemake.input, snakemake.output[0])
