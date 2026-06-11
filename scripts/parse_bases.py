import sys

sys.stderr = open(snakemake.log[0], "w")

input_file = snakemake.input.mpileup
output_file = snakemake.output.base_stat


with open(input_file) as f, open(output_file, 'w') as o:
    o.write("CHROM\tPOS\tREF\tA\tT\tC\tG\tTOTAL\n")
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 6:
            continue
        chrom, pos, ref, depth, bases, quals = parts[:6]

        counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        i = 0
        b_up = bases.upper()
        while i < len(b_up):
            b = b_up[i]
            if b in '.,':
                counts[ref.upper()] = counts.get(ref.upper(), 0) + 1
            elif b in 'ACGT':
                counts[b] += 1
            elif b == '^':
                i += 1          # 跳过后面的 mapping quality 字符
            elif b in '+-':     # indel，跳过插入/缺失的碱基
                i += 1
                num_str = ''
                while i < len(b_up) and b_up[i].isdigit():
                    num_str += b_up[i]
                    i += 1
                i += int(num_str) - 1
            # '$' / '*' / '<' / '>' 直接跳过
            i += 1

        total = sum(counts.values())
        o.write(f"{chrom}\t{pos}\t{ref.upper()}\t"
                f"{counts['A']}\t{counts['T']}\t{counts['C']}\t{counts['G']}\t{total}\n")
