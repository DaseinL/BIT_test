#!/usr/bin/env python3
import subprocess, os

def main():
    input_fasta = "10^4.fasta"
    blast_output = "blast_output.tsv"

    # 1. 构建数据库
    subprocess.run([
        "makeblastdb", "-in", input_fasta, "-dbtype", "nucl", "-out", input_fasta
    ], check=True)
    print("makeblastdb complete")

    # 2. 自比对
    subprocess.run([
        "blastn", "-query", input_fasta, "-db", input_fasta, "-out", blast_output,
        "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
        "-evalue", "1e-5", "-dust", "no", "-soft_masking", "false", "-num_threads", "4"
    ], check=True)
    print("blastn result:", blast_output)

    # 3. 初步过滤
    rows = []
    with open(blast_output, "r") as f:
        for line in f:
            c = line.strip().split('\t')
            if len(c) < 10:
                continue
            pident = float(c[2])
            length = int(c[3])
            qs = int(c[4])
            qe = int(c[5])
            ss = int(c[6])
            se = int(c[7])
            # 基础过滤
            if length < 900:
                continue
            if pident < 70.0:
                continue
            if qs == ss and qe == se:
                continue
            rows.append({
                "qstart": qs,
                "qend": qe,
                "sstart": ss,
                "send": se,
                "pident": pident,
                "length": length
            })

    # 排序：按 query 坐标排序
    rows.sort(key=lambda x: (x["qstart"], x["qend"]))

    groups = []
    used_regions = []  # 用来记录已使用的 query 和 subject 区域，避免重复分组

    # 可调阈值
    MAX_GAP = 50  # query 坐标最大合并距离
    MAX_LEN_DIFF = 1000  # 长度差的阈值（示例）
    MAX_PID_DIFF = 5.0  # 相似度差异阈值（示例）
    
    def is_overlapping(region1, region2):
        """判断两个区域是否有重叠"""
        return not (region1[1] < region2[0] or region1[0] > region2[1])

    for h in rows:
        qs, qe = h["qstart"], h["qend"]
        ss, se = h["sstart"], h["send"]
        pident, hlen = h["pident"], h["length"]

        # 构造当前比对的唯一标识（query region 和 subject region）
        query_key = (qs, qe)
        subject_key = (ss, se)

        # 检查当前比对的 query 和 subject 是否已经被使用过
        is_duplicate = False
        for grp in groups:
            # 判断当前的 query region 是否与现有分组的 query region 有重叠
            if is_overlapping((qs, qe), (grp["qstart"], grp["qend"])):
                # 如果有重叠区域，合并区域
                grp["qstart"] = min(grp["qstart"], qs)
                grp["qend"] = max(grp["qend"], qe)
                grp["subject_hits"].append({
                    "sstart": ss,
                    "send": se,
                    "pident": pident,
                    "length": hlen
                })
                grp["sum_len"] += hlen
                grp["sum_pid"] += pident
                grp["count"] += 1
                is_duplicate = True
                break

        if not is_duplicate:
            # 如果没有重复，创建一个新的分组
            groups.append({
                "qstart": qs,
                "qend": qe,
                "subject_hits": [{
                    "sstart": ss,
                    "send": se,
                    "pident": pident,
                    "length": hlen
                }],
                "sum_len": hlen,
                "sum_pid": pident,
                "count": 1
            })

    # 输出结果
    with open(input_fasta + "_grouped_repeats.txt", "w") as f:
        for i, grp in enumerate(groups, 1):
            f.write(f"Group {i}: Query region [{grp['qstart']}-{grp['qend']}]\n")
            for sh in grp["subject_hits"]:
                f.write(f"  - Subject region [{sh['sstart']}-{sh['send']}], Identity={sh['pident']:.2f}%, Length={sh['length']}\n")
            f.write("\n")

    print("group result directory: grouped_repeats.txt")


if __name__ == "__main__":
    main()
