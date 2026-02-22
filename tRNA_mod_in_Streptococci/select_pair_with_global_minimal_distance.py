#!/usr/bin/env python
# usage: python assign_min_sum.py input.tsv output.tsv
import sys
import math
from collections import defaultdict

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def main(infile, outfile):
    try:
        from scipy.optimize import linear_sum_assignment
    except Exception as e:
        print("This script requires scipy. Install it via `pip install scipy`.")
        sys.exit(1)

    # read file：record line of minimal value for each (key1, key2)
    # when (key1,key2) re-occurring, keep the line in which the value is smaller and occur first.
    pair_best_value = {}
    pair_first_line  = {}
    key1_order = {}   # record the order of key1, so the output is sorted and reproducible.
    order_counter = 0

    header = None
    with open(infile, "r", encoding="utf-8") as f:
        for idx, raw in enumerate(f):
            line = raw.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue

            # 首行若第三列不是数值，视为表头
            if idx == 0 and not is_float(cols[2]):
                header = line
                continue

            k1, k2 = cols[0], cols[1]
            vstr = cols[2]
            if not is_float(vstr):
                continue
            v = float(vstr)

            if k1 not in key1_order:
                key1_order[k1] = order_counter
                order_counter += 1

            key = (k1, k2)
            if key not in pair_best_value or v < pair_best_value[key]:
                pair_best_value[key] = v
                pair_first_line[key] = line  # 记录用于输出的整行（5列）

    # 构建 key1、key2 索引
    key1_list = sorted(key1_order.keys(), key=lambda k: key1_order[k])  # 按出现顺序
    key2_set = set(k2 for (_, k2) in pair_best_value.keys())
    key2_list = sorted(key2_set)  # 顺序任意即可

    nA = len(key1_list)
    nB = len(key2_list)

    if nA == 0:
        # 没有可用数据
        with open(outfile, "w", encoding="utf-8") as out:
            if header is not None:
                out.write(header + "\n")
        return

    if nB < nA:
        print(f"Impossible: |B| ({nB}) < |A| ({nA}) while requiring unique key2.")
        sys.exit(1)

    # 映射到索引
    a2i = {a:i for i,a in enumerate(key1_list)}
    b2j = {b:j for j,b in enumerate(key2_list)}

    # 构建代价矩阵（缺失边用大M）
    M = 1e15
    cost = [[M]*nB for _ in range(nA)]

    # 标记每个 key1 是否至少有一个可选 key2
    has_edge = [False]*nA

    for (k1, k2), v in pair_best_value.items():
        i = a2i[k1]
        j = b2j[k2]
        cost[i][j] = min(cost[i][j], v)
        has_edge[i] = True

    # 检查是否存在无可选 key2 的 key1
    for i, ok in enumerate(has_edge):
        if not ok:
            print(f"Impossible: key1 '{key1_list[i]}' has no candidate key2.")
            sys.exit(1)

    # 线性指派：行=key1，列=key2
    row_ind, col_ind = linear_sum_assignment(cost)

    # 过滤掉被迫选到 M 的不合法匹配（如果出现，说明无解）
    selected_pairs = []
    total_cost = 0.0
    for i, j in zip(row_ind, col_ind):
        c = cost[i][j]
        if c >= M/2:
            print("No feasible perfect matching: some key1 must take a forbidden key2.")
            sys.exit(1)
        a = key1_list[i]
        b = key2_list[j]
        selected_pairs.append((a, b))
        total_cost += c

    # 为输出找到对应的原始整行（用之前记录的 pair_first_line）
    # 注意：pair_first_line 里保存的是该 (a,b) 最小 value 的最早出现整行
    selected_lines = []
    # 为了输出稳定，按 key1 首次出现顺序输出
    selected_pairs.sort(key=lambda p: key1_order[p[0]])

    for a, b in selected_pairs:
        line = pair_first_line.get((a, b))
        if line is None:
            # 理论上不会发生，因为我们只对存在的边进行匹配
            line = f"{a}\t{b}\t{pair_best_value[(a,b)]}"
        selected_lines.append(line)

    with open(outfile, "w", encoding="utf-8") as out:
        if header is not None:
            out.write(header + "\n")
        for line in selected_lines:
            out.write(line + "\n")

    # 可选：在 stderr 打印总和（或写到文件尾部你也可以改）
    print(f"Selected {len(selected_lines)} lines. Min total of column 3 = {total_cost:.6f}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python assign_min_sum.py input.tsv output.tsv")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
