# -*- coding:utf-8 -*-
# @Time    :2023/3/23 18:14
# 2023-11-11 添加数据质控表：原始数据，umi标记数据，合并PCR产物数据统计
import os
import gzip
import time
import yaml
import datetime
import subprocess
import glob
import log
import shutil
import random  # 如果后面你想自己实现采样也可以用，这里先留着

logger = log.createCustomLogger('root')
def dep(file1, file2, config,out_dir):
    outfiles_r1 = {}
    outfiles_r2 = {}
    samples_name = {}
    for sample, value in config["samples"].items():
        if sample != "control":
            samples_name["{}".format(config["samples"]["{}".format(sample)]["barcode1"])] = sample
            outfiles_r1[sample] = open(os.path.join(out_dir, '%s.r1.tempumitagged.fastq' % sample), 'a')
            outfiles_r2[sample] = open(os.path.join(out_dir, '%s.r2.tempumitagged.fastq' % sample), 'a')
    
    #linker实验上用到的接头的前10个碱基，用来判断在哪一条read上出现linker，从而知道UMI位置是在R1还是在R2
    linker = "CTCCCTCGCC"  # The first 10 bases of linker used

    all_reads=0
    all_umi_reads=0

    for file31, file32 in zip(file1, file2):
        f31 = gzip.open(file31, 'rb')
        f32 = gzip.open(file32, 'rb')
        s1_1 = f31.readline().decode('utf-8')
        s1_2 = f31.readline().decode('utf-8')
        s1_3 = f31.readline().decode('utf-8')
        s1_4 = f31.readline().decode('utf-8')

        s2_1 = f32.readline().decode('utf-8')
        s2_2 = f32.readline().decode('utf-8')
        s2_3 = f32.readline().decode('utf-8')
        s2_4 = f32.readline().decode('utf-8')
        # .replace("\n", "")[::-1] + "\n"
        while s1_1:
            all_reads+=1
            # 情况1：linker是在R1上，barcode也在R1上。
            # s1_2[34:49]是查找R1序列的第35-49个碱基区间，在这里找linker
            # s1_2[:8]是查找R1序列的前8个碱基，在这里找barcode
            if linker in s1_2[34:49] and s1_2[:8] in samples_name.keys():
                all_umi_reads+=1
                outfiles_r1['{}'.format(samples_name["{}".format(s1_2[:8])])].write(
                    # s1_2[26:34]：R1 的第 27~34 bp（8bp）+ s1_2[46:52]：R1 的第 47~52 bp（6bp）+ s2_2[:6]：R2 的前 6bp =20bp的UMI
                    # 此时，ID行为：原ID（去掉换行，最后有个空格） + UMI + 换行 . 原来：@NB552061:123:HGYYGBGXH:1:11101:10000:1000 1:N:0:1\n 变成：@NB552061:... 1:N:0:1 ACGTACGT_TTGGCC_AAAAAA\n
                    # 序列行：s1_2[46:]把 R1 序列从第 47 个碱基开始保留下来，前面 46bp（里面包含 barcode、部分 UMI、linker 等）都丢掉。
                    # R2的ID同R1加上UMI，但序列保留完整
                    s1_1.replace("\n", " ") + s1_2[26:34] + "_" + s1_2[46:52] + "_" + s2_2[:6] + "\n" + s1_2[
                                                                                                        46:] + s1_3 + s1_4[
                                                                                                                      46:])
                outfiles_r2['{}'.format(samples_name["{}".format(s1_2[:8])])].write(
                    s2_1.replace("\n", " ") + s1_2[26:34] + "_" + s1_2[46:52] + "_" + s2_2[
                                                                                      :6] + "\n" + s2_2 + s2_3 + s2_4)
            # 情况2：R1 上没 linker，R2 上有 linker，barcode 在 R2
            elif linker not in s1_2[34:49] and linker in s2_2[34:49] and s2_2[:8] in samples_name.keys():
                all_umi_reads+=1

                outfiles_r1['{}'.format(samples_name["{}".format(s2_2[:8])])].write(
                    s2_1.replace("\n", " ") + s2_2[26:34] + "_" + s2_2[46:52] + "_" + s1_2[:6] + "\n" + s2_2[
                                                                                                        46:] + s2_3 + s2_4[
                                                                                                                      46:])
                outfiles_r2['{}'.format(samples_name["{}".format(s2_2[:8])])].write(
                    s1_1.replace("\n", " ") + s2_2[26:34] + "_" + s2_2[46:52] + "_" + s1_2[
                                                                                      :6] + "\n" + s1_2 + s1_3 + s1_4)
            # 情况3：如果两个条件都不满足（比如没找到 linker、或者 barcode 不是任何样本的），这一对 read 就被丢弃，不输出，也不会计入 all_umi_reads
            s1_1 = f31.readline().decode('utf-8')
            s1_2 = f31.readline().decode('utf-8')
            s1_3 = f31.readline().decode('utf-8')
            s1_4 = f31.readline().decode('utf-8')

            s2_1 = f32.readline().decode('utf-8')
            s2_2 = f32.readline().decode('utf-8')
            s2_3 = f32.readline().decode('utf-8')
            s2_4 = f32.readline().decode('utf-8')
        f31.close()
        f32.close()
    for sample, value in config["samples"].items():
        if sample != "control":
            samples_name["{}".format(config["samples"]["{}".format(sample)]["barcode1"])] = sample
            outfiles_r1[sample].close()
            outfiles_r2[sample].close()

    # 对每个样本的输出fastq按UMI排序
    for sample, value in config["samples"].items():
        if sample != "control":
            r1_umitagged_unsorted_file = os.path.join(out_dir, '%s.r1.tempumitagged.fastq' % sample)
            r2_umitagged_unsorted_file = os.path.join(out_dir, '%s.r2.tempumitagged.fastq' % sample)
            read1_out = os.path.join(out_dir, '%s.r1.umitagged.fastq' % sample)
            read2_out = os.path.join(out_dir, '%s.r2.umitagged.fastq' % sample)
            # 其中的paste- - - -是将fastq的四行合并，便于sort对整条read操作；tr "\t" "\n"再拆回原来的fastq格式内容
            # sort -k3,3 -k1,1:按照第3列优先，再按第1列排序（排序的目的是将同一个UMI的reads聚在一起）
            cmd = 'cat ' + r1_umitagged_unsorted_file + ' | paste - - - - | sort -k3,3 -k1,1 | tr "\t" "\n" >' + read1_out
            subprocess.check_call(cmd, shell=True, env=os.environ.copy())
            cmd = 'cat ' + r2_umitagged_unsorted_file + ' | paste - - - - | sort -k3,3 -k1,1 | tr "\t" "\n" >' + read2_out
            subprocess.check_call(cmd, shell=True, env=os.environ.copy())
            # 删除未进行排序的文件
            os.remove(r1_umitagged_unsorted_file)
            os.remove(r2_umitagged_unsorted_file)

    return all_reads,all_umi_reads

def count_fastq_reads(fastq_path):
    """
    统计 FASTQ 文件中的 read 数量。每 4 行是一条 read。
    """
    line_count = 0
    with open(fastq_path, 'r') as f:
        for _ in f:
            line_count += 1
    if line_count == 0:
        return 0
    return line_count // 4

def normalize_umitagged(out_dir="./umitagged", target_count=100000, seed=100):
    """
    对 out_dir 下的 *.r1.umitagged.fastq / *.r2.umitagged.fastq 做下采样，
    生成 out_dir_normalized 目录，然后把原目录改名为 out_dir_original，
    再把 normalized 目录改名为 out_dir。
    同时统计每个样本下采样前后的 reads 数，写到 reads_normalization.txt。
    """
    logger.info(f"Start normalization / downsampling in {out_dir}, target_count={target_count}")
    normalized_dir = out_dir + "_normalized"
    original_dir = out_dir + "_original"

    if not os.path.exists(normalized_dir):
        os.mkdir(normalized_dir)

    # 遍历所有 r1.umitagged.fastq
    r1_files = glob.glob(os.path.join(out_dir, "*r1.umitagged.fastq"))
    if not r1_files:
        logger.warning(f"No *r1.umitagged.fastq found in {out_dir}, skip normalization.")
        return

    # 用来记录每个样本的 reads 统计信息
    # 格式：(sample_name, r1_before, r2_before, r1_after, r2_after)
    stats = []

    for r1_file in r1_files:
        sample_name = os.path.basename(r1_file).replace(".r1.umitagged.fastq", "")
        r2_file = os.path.join(out_dir, f"{sample_name}.r2.umitagged.fastq")
        if not os.path.exists(r2_file):
            logger.warning(f"R2 file not found for sample {sample_name}, skip this sample.")
            continue

        # 统计 read 数（采样前）
        r1_count = count_fastq_reads(r1_file)
        r2_count = count_fastq_reads(r2_file)

        if r1_count == 0 or r2_count == 0:
            logger.warning(f"One of the files for {sample_name} has 0 reads. Skipping this sample.")
            continue

        if r1_count != r2_count:
            logger.warning(f"R1/R2 read count mismatch for {sample_name}: R1={r1_count}, R2={r2_count}")

        logger.info(f"Sample {sample_name}: R1={r1_count}, R2={r2_count}")

        # 决定要保留多少条
        if r1_count < target_count:
            # 序列数不足 target 就全保留（数量不变，只是随机重排）
            keep_n = r1_count
        else:
            keep_n = target_count

        out_r1 = os.path.join(normalized_dir, f"{sample_name}.r1.umitagged.fastq")
        out_r2 = os.path.join(normalized_dir, f"{sample_name}.r2.umitagged.fastq")

        # 使用 seqtk 做随机抽样（保持 FASTQ 格式）
        cmd1 = f"seqtk sample -s{seed} {r1_file} {keep_n} > {out_r1}"
        cmd2 = f"seqtk sample -s{seed} {r2_file} {keep_n} > {out_r2}"

        logger.info(f"Downsampling {sample_name} R1 to {keep_n}: {cmd1}")
        subprocess.check_call(cmd1, shell=True, env=os.environ.copy())

        logger.info(f"Downsampling {sample_name} R2 to {keep_n}: {cmd2}")
        subprocess.check_call(cmd2, shell=True, env=os.environ.copy())

        logger.info(f"Normalized {sample_name} to {keep_n} reads")

        # 记录统计：采样前 / 采样后
        stats.append((sample_name, r1_count, r2_count, keep_n, keep_n))

    # 目录重命名：umitagged -> umitagged_original, umitagged_normalized -> umitagged
    if os.path.exists(original_dir):
        logger.warning(f"{original_dir} already exists, removing it.")
        shutil.rmtree(original_dir)

    logger.info(f"Renaming {out_dir} -> {original_dir}")
    os.rename(out_dir, original_dir)

    logger.info(f"Renaming {normalized_dir} -> {out_dir}")
    os.rename(normalized_dir, out_dir)

    logger.info("Normalization and renaming completed.")

    # 写出采样前后的 reads 统计信息
    # 注意：此时 out_dir 已经指向“新的归一化后的 umitagged 目录”
    stats_file = os.path.join(out_dir, "reads_normalization.txt")
    with open(stats_file, "w") as sf:
        sf.write("sample\tr1_before\tr2_before\tr1_after\tr2_after\n")
        for sample_name, r1_b, r2_b, r1_a, r2_a in stats:
            sf.write(f"{sample_name}\t{r1_b}\t{r2_b}\t{r1_a}\t{r2_a}\n")

    logger.info(f"Normalization stats written to {stats_file}")

# 发现连续相同的UMI的reads，选取质量总和最高的read作为代表，输出新read，ID变成@UMI_重复次数，表示这UMI对应多少条PCR dupliacte集合
def consolidate(file_umi):
    logger.info("loading {}".format(file_umi))

    pcr_con_reads=0

    if not os.path.exists("./consolidated"):
        os.mkdir("./consolidated")
    outf = file_umi.replace("umitagged.", "consolidated.").replace("umitagged", "consolidated")
    outfile = open(outf, 'a')
    with open(file_umi, 'r') as f:
        s1 = f.readline()
        s2 = f.readline()
        s3 = f.readline()
        s4 = f.readline()
        # front_umi_id：当前正在积累的这一组 UMI 的 ID；count就是这一组UMI内累计了多少条reads
        front_umi_id, count = "", 0

        while s1:
            cur_umi = s1.split(" ")[-1].replace("\n", "")
            # 第一次遇到read就初始化当前UMI组
            # front_seq当前组里目前最好的序列；front_q对应的质量行；all_q当前这条read的总质量分数
            # ord(x) - 33 ：从 ASCII 字符计算 Phred 质量值；求和表示“这条 read 的整体质量有多好”。
            if front_umi_id == "":
                front_seq, front_q, all_q = s2, s4, sum([(ord(x) - 33) for x in s4.replace("\n", "")])
                front_umi_id = s1.split(" ")[-1].replace("\n", "")
                count += 1
                # 读下一条
                s1 = f.readline()
                s2 = f.readline()
                s3 = f.readline()
                s4 = f.readline()
                continue
            # 遇到新的UMI，先将上一组的UMI的代表的read写出去
            if cur_umi != front_umi_id:
                pcr_con_reads+=1
                outfile.write("@" + front_umi_id + "_{}\n".format(count) + front_seq + s3 + front_q)
                front_umi_id = cur_umi
                front_seq, front_q, all_q, count = s2, s4, sum([(ord(x) - 33) for x in s4.replace("\n", "")]), 1
            # 如何还是同一个UMI，计算当前read的质量总和，如果大于前一个，将这条read替换成当前UMI组的最佳代表，count+1
            else:
                cur_q = sum([(ord(x) - 33) for x in s4.replace("\n", "")])
                if cur_q > all_q:
                    front_seq, front_q, all_q = s2, s4, cur_q
                count += 1
            s1 = f.readline()
            s2 = f.readline()
            s3 = f.readline()
            s4 = f.readline()
        f.close()
    outfile.close()
    logger.info("the result has been saved {}".format(outfile))

    return pcr_con_reads


def main(config,data1,data2):  # manifest_data
    logger.info(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    out_dir = "./umitagged"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    file_1 = data1
    file_2 = data2
    all_reads,all_umi_reads = dep(file_1, file_2, config,out_dir)
    # 2. 在 ./umitagged 基础上做下采样 / 归一化，并重命名目录
    normalize_umitagged(out_dir=out_dir, target_count=100000, seed=100)

    all_pcr_reads=0
    file_umi = glob.glob("./umitagged/*.fastq")
    for fileumi in file_umi:
        all_pcr_reads+=consolidate(file_umi=fileumi)
    all_pcr_reads = all_pcr_reads/2
    umi_by_all = float(all_umi_reads)/all_reads
    pcr_by_all = float(all_pcr_reads)/all_reads

    with open("reads_statistics.txt","w") as rf:
        rf.write("all reads:{} (1) \nall umi reads:{} ({}) \nall pcr reads:{} ({})".format(all_reads, all_umi_reads, umi_by_all, all_pcr_reads, pcr_by_all))
    logger.info(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
