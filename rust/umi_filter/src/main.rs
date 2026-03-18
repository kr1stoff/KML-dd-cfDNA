// 导入标准库模块
use std::collections::{HashMap, HashSet}; // 用于存储 UMI 序列
use std::env; // 用于获取命令行参数
use std::fs::File; // 用于文件操作
use std::io::{BufRead, BufReader, BufWriter, Write}; // 用于文件读写

// 常量定义
const MAX_MISMATCH: usize = 1; // 允许的最大错配数
const UMI_LEN: usize = 6; // UMI 序列长度
const Q20_THRESHOLD: u8 = 20; // Q20 质量值阈值
const Q30_THRESHOLD: u8 = 30; // Q30 质量值阈值
const BUFFER_SIZE: usize = 1024 * 1024; // 设置 1MB 读取缓冲区

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

/// 统计质量值字符串中 >= 指定阈值的碱基数
///
/// # 参数
/// * `qual_str` - 质量值字符串(Phred33 编码)
/// * `threshold` - 质量值阈值
///
/// # 返回值
/// 满足条件的碱基数
fn count_quality(qual_str: &str, threshold: u8) -> usize {
    qual_str
        .chars()
        .filter(|&c| {
            let quality = c as u8 - 33; // Phred33 编码转换
            quality >= threshold
        })
        .count()
}

/// 从文件读取 UMI 序列列表
///
/// # 参数
/// * `umi_file_path` - UMI 文件路径
///
/// # 返回值
/// 成功时返回 UMI 序列列表，失败时返回 io::Error
fn read_umi_list(umi_file_path: &str) -> std::io::Result<Vec<String>> {
    let file = File::open(umi_file_path)?;
    let reader = BufReader::new(file);

    let mut umi_list = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if !trimmed.is_empty() {
            umi_list.push(trimmed.to_string());
        }
    }

    Ok(umi_list)
}

/// 生成允许的 UMI 序列集合
///
/// 该函数会：
/// 1. 生成所有 UMI 序列的组合(每个 UMI 与其他所有 UMI 组合)
/// 2. 为每个组合生成允许 1 个错配的变体
/// 3. 返回一个 HashMap，键为允许的序列，值为对应的原始序列
///
/// # 参数
/// * `umi_list` - 原始 UMI 序列列表
///
/// # 返回值
/// 包含所有允许的 UMI 序列及其对应原始序列的 HashMap
fn generate_allowed(umi_list: &[String]) -> HashMap<String, String> {
    // 定义核苷酸可能的碱基
    let nucs = ['A', 'C', 'G', 'T'];

    // 预分配 whitelist 容量，避免频繁扩容
    let whitelist_capacity = umi_list.len() * umi_list.len();
    let mut whitelist = HashSet::with_capacity(whitelist_capacity);

    // 生成所有可能的 UMI 组合(u1 + u2)
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
                        // 创建新的序列(替换一个碱基)
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
    println!("允许序列总数(含1 mismatch邻居): {}", allowed.len());

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
/// 6. 统计并输出过滤结果到文件
///
/// # 参数
/// * `r1_path` - R1 输入文件路径
/// * `r2_path` - R2 输入文件路径
/// * `r1_out_path` - R1 输出文件路径
/// * `r2_out_path` - R2 输出文件路径
/// * `allowed` - 允许的 UMI 序列集合
/// * `stats_out_path` - 统计信息输出文件路径
///
/// # 返回值
/// 成功时返回 Ok(())，失败时返回 io::Error
fn filter_fastq(
    r1_path: &str,
    r2_path: &str,
    r1_out_path: &str,
    r2_out_path: &str,
    allowed: &HashMap<String, String>,
    stats_out_path: &str,
) -> std::io::Result<()> {
    // 打开未压缩的输入文件
    let r1_file = File::open(r1_path)?;
    let r2_file = File::open(r2_path)?;

    // 创建缓冲读取器
    let r1_reader = BufReader::with_capacity(BUFFER_SIZE, r1_file);
    let r2_reader = BufReader::with_capacity(BUFFER_SIZE, r2_file);

    // 创建未压缩的输出文件
    let out1 = File::create(r1_out_path)?;
    let out2 = File::create(r2_out_path)?;

    // 创建缓冲写入器
    let mut out1 = BufWriter::with_capacity(BUFFER_SIZE, out1);
    let mut out2 = BufWriter::with_capacity(BUFFER_SIZE, out2);

    // 创建行迭代器
    let mut r1_lines = r1_reader.lines();
    let mut r2_lines = r2_reader.lines();

    // 统计变量
    let mut total_reads = 0usize; // 总读取数
    let mut kept_reads = 0usize; // 保留的读取数
    let mut total_bases = 0usize; // 总碱基数
    let mut q20_count = 0usize; // Q20 碱基数
    let mut q30_count = 0usize; // Q30 碱基数

    // 预分配 UMI 组合字符串，避免频繁分配
    let mut combined = String::with_capacity(UMI_LEN * 2);

    // 循环处理每对读取
    loop {
        // 读取 R1 的第一行(标题行)
        let r1 = match r1_lines.next() {
            Some(v) => v?, // 解析结果，处理可能的错误
            None => break, // 文件结束，退出循环
        };

        // 读取 R1 的剩余三行(序列、+、质量值)
        let r1_seq = r1_lines.next().unwrap()?;
        let r1_plus = r1_lines.next().unwrap()?;
        let r1_qual = r1_lines.next().unwrap()?;

        // 读取 R2 的四行
        let r2 = r2_lines.next().unwrap()?;
        let r2_seq = r2_lines.next().unwrap()?;
        let r2_plus = r2_lines.next().unwrap()?;
        let r2_qual = r2_lines.next().unwrap()?;

        // 统计质量值
        q20_count += count_quality(&r1_qual, Q20_THRESHOLD);
        q20_count += count_quality(&r2_qual, Q20_THRESHOLD);
        q30_count += count_quality(&r1_qual, Q30_THRESHOLD);
        q30_count += count_quality(&r2_qual, Q30_THRESHOLD);

        // 统计读取数和碱基数
        total_reads += 1;
        total_bases += r1_seq.len() + r2_seq.len();

        // 提取 UMI 序列并组合
        combined.clear();
        combined.push_str(extract_umi(&r1_seq));
        combined.push_str(extract_umi(&r2_seq));

        // 检查是否在允许列表中
        if allowed.contains_key(&combined) {
            // 写入符合条件的读取到输出文件(切掉 UMI 部分)
            writeln!(out1, "{}", r1)?;
            writeln!(out1, "{}", &r1_seq[UMI_LEN..])?; // 切掉前 6bp UMI
            writeln!(out1, "{}", r1_plus)?;
            writeln!(out1, "{}", &r1_qual[UMI_LEN..])?; // 切掉前 6bp 质量值

            writeln!(out2, "{}", r2)?;
            writeln!(out2, "{}", &r2_seq[UMI_LEN..])?; // 切掉前 6bp UMI
            writeln!(out2, "{}", r2_plus)?;
            writeln!(out2, "{}", &r2_qual[UMI_LEN..])?; // 切掉前 6bp 质量值

            // 增加保留计数
            kept_reads += 1;
        }
    }

    // 计算 Q20 和 Q30 比例
    let q20_ratio = if total_bases > 0 {
        q20_count as f64 / total_bases as f64
    } else {
        0.0
    };

    let q30_ratio = if total_bases > 0 {
        q30_count as f64 / total_bases as f64
    } else {
        0.0
    };

    // 输出统计信息到文件
    let mut stats_file = File::create(stats_out_path)?;
    writeln!(stats_file, "total_reads: {}", total_reads)?;
    writeln!(stats_file, "kept_reads: {}", kept_reads)?;
    writeln!(stats_file, "total_bases: {}", total_bases)?;
    writeln!(stats_file, "q20: {:.4}", q20_ratio)?;
    writeln!(stats_file, "q30: {:.4}", q30_ratio)?;

    // 打印过滤结果
    println!(
        "总 reads: {}, 保留 reads: {}, 过滤掉: {}",
        total_reads,
        kept_reads,
        total_reads - kept_reads
    );
    println!("Q20 比例: {:.4}, Q30 比例: {:.4}", q20_ratio, q30_ratio);

    Ok(())
}

/// 程序主函数
///
/// 该函数会：
/// 1. 解析命令行参数
/// 2. 从文件读取 UMI 序列列表
/// 3. 生成允许的 UMI 序列集合
/// 4. 调用过滤函数处理 FASTQ 文件
///
/// # 返回值
/// 成功时返回 Ok(())，失败时返回 io::Error
fn main() -> std::io::Result<()> {
    // 获取命令行参数
    let args: Vec<String> = env::args().collect();

    // 检查参数数量是否正确
    if args.len() != 7 {
        eprintln!("Usage: umi_filter R1_in.fastq R2_in.fastq R1_out.fastq R2_out.fastq umi_file stats_outfile");
        std::process::exit(1);
    }

    // 从文件读取 UMI 序列列表
    let umi_list = read_umi_list(&args[5])?;
    println!("从文件读取到 {} 个 UMI 序列", umi_list.len());

    // 生成允许的 UMI 序列集合
    let allowed = generate_allowed(&umi_list);

    // 调用过滤函数处理 FASTQ 文件
    filter_fastq(&args[1], &args[2], &args[3], &args[4], &allowed, &args[6])?;

    Ok(())
}
