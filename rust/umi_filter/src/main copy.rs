// 导入标准库模块
use std::collections::{HashMap, HashSet}; // 用于存储 UMI 序列
use std::env; // 用于获取命令行参数
use std::fs::File; // 用于文件操作
use std::io::{BufRead, BufReader, BufWriter, Write}; // 用于文件读写

// 常量定义
const MAX_MISMATCH: usize = 1; // 允许的最大错配数
const UMI_LEN: usize = 6; // UMI 序列长度

/// 从序列中提取 UMI 部分
///
/// # 参数
/// * `seq` - 原始序列字符串
///
/// # 返回值
/// 提取的 UMI 序列
fn extract_umi(seq: &str) -> &str {
    &seq[..UMI_LEN]
}

/// 生成允许的 UMI 序列集合
///
/// 该函数会：
/// 1. 生成所有 UMI 序列的组合（每个 UMI 与其他所有 UMI 组合）
/// 2. 为每个组合生成允许 1 个错配的变体
/// 3. 返回一个 HashMap，键为允许的序列，值为对应的原始序列
///
/// # 参数
/// * `umi_list` - 原始 UMI 序列列表
///
/// # 返回值
/// 包含所有允许的 UMI 序列及其对应原始序列的 HashMap
fn generate_allowed(umi_list: &[&str]) -> HashMap<String, String> {
    // 定义核苷酸可能的碱基
    let nucs = ['A', 'C', 'G', 'T'];

    // 预分配 whitelist 容量，避免频繁扩容
    let whitelist_capacity = umi_list.len() * umi_list.len();
    let mut whitelist = HashSet::with_capacity(whitelist_capacity);

    // 生成所有可能的 UMI 组合（u1 + u2）
    for u1 in umi_list {
        for u2 in umi_list {
            whitelist.insert(format!("{}{}", u1, u2));
        }
    }

    // 预分配 allowed 容量，估计为原始序列数加上可能的错配数
    let estimated_capacity = whitelist.len() * (1 + UMI_LEN * 3); // 原始 + 每个位置3个可能的错配
    let mut allowed = HashMap::with_capacity(estimated_capacity);

    // 为每个 UMI 组合及其变体添加到允许列表
    for w in &whitelist {
        // 添加原始序列
        allowed.insert(w.clone(), w.clone());

        // 如果允许 1 个错配
        if MAX_MISMATCH >= 1 {
            // 将序列转换为字符向量
            let chars: Vec<char> = w.chars().collect();

            // 遍历每个位置
            for i in 0..chars.len() {
                // 尝试替换为其他碱基
                for n in &nucs {
                    if *n != chars[i] {
                        // 创建新的序列（替换一个碱基）
                        let mut new_chars = chars.clone();
                        new_chars[i] = *n;
                        let neighbor: String = new_chars.iter().collect();
                        // 添加到允许列表，值为原始序列
                        allowed.insert(neighbor, w.clone());
                    }
                }
            }
        }
    }

    // 打印允许的序列总数
    println!("允许序列总数（含1 mismatch邻居）: {}", allowed.len());

    allowed
}

/// 过滤 FASTQ 文件，只保留含有允许 UMI 序列的读取
///
/// 该函数会：
/// 1. 打开输入的未压缩 FASTQ 文件
/// 2. 创建输出的未压缩 FASTQ 文件
/// 3. 逐行读取输入文件，提取 UMI 序列
/// 4. 检查 UMI 序列是否在允许列表中
/// 5. 将符合条件的读取写入输出文件
/// 6. 统计并打印过滤结果
///
/// # 参数
/// * `r1_path` - R1 输入文件路径
/// * `r2_path` - R2 输入文件路径
/// * `r1_out_path` - R1 输出文件路径
/// * `r2_out_path` - R2 输出文件路径
/// * `allowed` - 允许的 UMI 序列集合
///
/// # 返回值
/// 成功时返回 Ok(())，失败时返回 io::Error
fn filter_fastq(
    r1_path: &str,
    r2_path: &str,
    r1_out_path: &str,
    r2_out_path: &str,
    allowed: &HashMap<String, String>,
) -> std::io::Result<()> {
    // 打开未压缩的输入文件
    let r1_file = File::open(r1_path)?;
    let r2_file = File::open(r2_path)?;

    // 创建缓冲读取器
    let r1_reader = BufReader::new(r1_file);
    let r2_reader = BufReader::new(r2_file);

    // 创建未压缩的输出文件
    let out1 = File::create(r1_out_path)?;
    let out2 = File::create(r2_out_path)?;

    // 创建缓冲写入器
    let mut out1 = BufWriter::new(out1);
    let mut out2 = BufWriter::new(out2);

    // 创建行迭代器
    let mut r1_lines = r1_reader.lines();
    let mut r2_lines = r2_reader.lines();

    // 统计变量
    let mut total = 0usize; // 总读取数
    let mut kept = 0usize; // 保留的读取数

    // 预分配 UMI 组合字符串，避免频繁分配
    let mut combined = String::with_capacity(UMI_LEN * 2);

    // 循环处理每对读取
    loop {
        // 读取 R1 的第一行（标题行）
        let r1 = match r1_lines.next() {
            Some(v) => v?, // 解析结果，处理可能的错误
            None => break, // 文件结束，退出循环
        };

        // 读取 R1 的剩余三行（序列、+、质量值）
        let r1_seq = r1_lines.next().unwrap()?;
        let r1_plus = r1_lines.next().unwrap()?;
        let r1_qual = r1_lines.next().unwrap()?;

        // 读取 R2 的四行
        let r2 = r2_lines.next().unwrap()?;
        let r2_seq = r2_lines.next().unwrap()?;
        let r2_plus = r2_lines.next().unwrap()?;
        let r2_qual = r2_lines.next().unwrap()?;

        // 增加总读取数
        total += 1;

        // 提取 UMI 序列并组合
        combined.clear();
        combined.push_str(extract_umi(&r1_seq));
        combined.push_str(extract_umi(&r2_seq));

        // 检查是否在允许列表中
        if allowed.contains_key(&combined) {
            // 写入符合条件的读取到输出文件
            writeln!(out1, "{}", r1)?;
            writeln!(out1, "{}", r1_seq)?;
            writeln!(out1, "{}", r1_plus)?;
            writeln!(out1, "{}", r1_qual)?;

            writeln!(out2, "{}", r2)?;
            writeln!(out2, "{}", r2_seq)?;
            writeln!(out2, "{}", r2_plus)?;
            writeln!(out2, "{}", r2_qual)?;

            // 增加保留计数
            kept += 1;
        }
    }

    // 打印过滤结果
    println!(
        "总 reads: {}, 保留 reads: {}, 过滤掉: {}",
        total,
        kept,
        total - kept
    );

    Ok(())
}

/// 程序主函数
///
/// 该函数会：
/// 1. 解析命令行参数
/// 2. 定义 UMI 序列列表
/// 3. 生成允许的 UMI 序列集合
/// 4. 调用过滤函数处理 FASTQ 文件
///
/// # 返回值
/// 成功时返回 Ok(())，失败时返回 io::Error
fn main() -> std::io::Result<()> {
    // 获取命令行参数
    let args: Vec<String> = env::args().collect();

    // 检查参数数量是否正确
    if args.len() != 5 {
        eprintln!("Usage: umi_filter R1_in.fastq R2_in.fastq R1_out.fastq R2_out.fastq");
        std::process::exit(1);
    }

    // 定义原始 UMI 序列列表
    let umi_list = vec![
        "TCGGTA", "GACGTT", "AAGTTG", "AAGGCA", "AAGCGT", "AACTGA", "AACCAG", "ATAAGC", "ATACCT",
        "ATTGAG", "ATGCAA", "ATCATG", "ATCTAC", "AGATGT", "AGACTC", "AGTACA", "AGTGGC", "ACAACG",
        "ACAGAT", "TAATCG", "TATGAC", "TATCCA", "TACTTC", "TTACGA", "TTGTCA", "TTCACT", "TTCTGG",
        "TGAGCA", "TGTATC", "TGTTAG", "TGTCGT", "TGCCAC", "TCAGGC", "TCTTGA", "TCTGCT", "TCGACC",
        "TCGCAG", "TCCTAT", "GAATAC", "GATGGA", "GACAGC", "GTAGAA", "GTTAAC", "GTTCTA", "GTGAGT",
        "GGAACC", "GGATTA", "GGCAAG", "GCAATT", "GCACCA", "GCTTCG", "GCTCAT", "GCGTAA", "CAAGTA",
        "CAACAT", "CATTGT", "CATGCG", "CTAGGT", "CTGATA", "CTGGAC", "CTCCTC", "CGTCTG", "CCATAG",
        "CCGCTT",
    ];

    // 生成允许的 UMI 序列集合
    let allowed = generate_allowed(&umi_list);

    // 调用过滤函数处理 FASTQ 文件
    filter_fastq(&args[1], &args[2], &args[3], &args[4], &allowed)?;

    Ok(())
}
