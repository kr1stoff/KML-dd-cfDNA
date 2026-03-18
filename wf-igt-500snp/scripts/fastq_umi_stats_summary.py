import sys

sys.stderr = open(snakemake.log[0], "w")
