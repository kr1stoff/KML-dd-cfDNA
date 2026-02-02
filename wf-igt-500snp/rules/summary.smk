rule all_qc_summary:
    input:
        rules.fq_stats_summary.output,
        rules.bam_stats_summary.output,
    output:
        "upload/panel-qc-summary.tsv",
        "upload/panel-qc-summary.xlsx",
    benchmark:
        ".log/summary/all_qc_summary.bm"
    log:
        ".log/summary/all_qc_summary.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/all_qc_summary.py"
