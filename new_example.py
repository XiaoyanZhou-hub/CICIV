#!/usr/bin/env python
# encoding: utf-8

import sys

sys.path.append(r"C:\Users\62589\OneDrive\Desktop\Bib-jiaojie\PC-simple\pyCausalFS\CBD")

from MBs.pc_simple import pc_simple
from MBs.MMMB.MMPC import MMPC
from MBs.HITON.HITON_PC import HITON_PC
import tarfile
import pandas as pd
import numpy as np
import os
import time
from scipy import stats

# 40个目标基因（按TCGA因果排序）
TARGET_GENES = [
    "FOXA1", "AKT1", "NTRK3", "KEAP1", "GOLPH3", "IRS4", "EP300", "ESR1", "GATA3", "FBLN2",
    "IKZF3", "SALL4", "NOTCH1", "FADD", "CCND1", "MAP2K4", "BAP1", "ASPM", "PPM1D", "PBRM1",
    "PIK3CA", "TP53", "MAP3K1", "CTCF", "TBX3", "ETV6", "MAP3K13", "CDH1", "FLNA", "CASP8",
    "BRCA2", "SMARCD1", "BARD1", "RB1", "ZMYM3", "NCOR1", "CDKN1B", "ARID1B", "ARID1A", "ERBB2"
]

TCGA_ATE = {
    "FOXA1": 0.3005, "AKT1": 0.1715, "NTRK3": -0.1362, "KEAP1": 0.1302, "GOLPH3": 0.1217,
    "IRS4": -0.1210, "EP300": 0.1192, "ESR1": 0.1166, "GATA3": 0.1164, "FBLN2": -0.1139,
    "IKZF3": 0.1123, "SALL4": 0.1106, "NOTCH1": -0.1094, "FADD": 0.1053, "CCND1": 0.1037,
    "MAP2K4": 0.098, "BAP1": 0.0973, "ASPM": -0.0947, "PPM1D": -0.0896, "PBRM1": -0.0891,
    "PIK3CA": -0.0883, "TP53": -0.0854, "MAP3K1": 0.0815, "CTCF": 0.0808, "TBX3": -0.0798,
    "ETV6": -0.0794, "MAP3K13": 0.0793, "CDH1": 0.0756, "FLNA": -0.0692, "CASP8": 0.0659,
    "BRCA2": 0.0655, "SMARCD1": -0.0635, "BARD1": 0.0626, "RB1": -0.0623, "ZMYM3": -0.0573,
    "NCOR1": -0.0538, "CDKN1B": 0.0525, "ARID1B": -0.0512, "ARID1A": 0.0421, "ERBB2": 0.0333
}


def extract_metabric(tar_path, extract_dir="metabric_data"):
    if not os.path.exists(extract_dir):
        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extractall(extract_dir)
    return extract_dir


def load_metabric_data(data_dir):
    expr_file = os.path.join(data_dir, "brca_metabric", "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt")
    clinical_file = os.path.join(data_dir, "brca_metabric", "data_clinical_patient.txt")

    expr_df = pd.read_csv(expr_file, sep="\t", index_col=0, low_memory=False)
    if "Entrez_Gene_Id" in expr_df.columns:
        expr_df = expr_df.drop("Entrez_Gene_Id", axis=1)
    expr_df = expr_df.T
    expr_df = expr_df.loc[:, ~expr_df.columns.duplicated()]

    clinical_df = pd.read_csv(clinical_file, sep="\t", comment="#").set_index("PATIENT_ID")
    common = expr_df.index.intersection(clinical_df.index)
    expr_df = expr_df.loc[common]
    clinical_df = clinical_df.loc[common]

    survival = clinical_df["OS_STATUS"].apply(
        lambda x: 1 if "DECEASED" in str(x).upper() or ":1" in str(x) else 0
    )

    return expr_df, survival, clinical_df

# 离散化函数
def safe_discretize(arr, n_bins=3):

    result = np.zeros(len(arr), dtype=np.int32)

    # 转换为数值并处理NaN
    arr = pd.to_numeric(pd.Series(arr), errors='coerce')

    # 填充NaN为中位数
    if arr.isna().any():
        arr = arr.fillna(arr.median())

    if arr.isna().any():
        arr = arr.fillna(0)

    try:
        result = pd.qcut(arr.rank(method='first'), q=n_bins, labels=False).values.astype(np.int32)
    except:
        try:
            result = pd.cut(arr, bins=n_bins, labels=False).fillna(1).values.astype(np.int32)
        except:
            q1 = arr.quantile(0.33)
            q2 = arr.quantile(0.67)
            result = np.where(arr <= q1, 0, np.where(arr <= q2, 1, 2)).astype(np.int32)

    # 确保所有值都在[0, n_bins-1]范围内
    result = np.clip(result, 0, n_bins - 1).astype(np.int32)

    return result

# 对整个DataFrame进行离散化
def discretize_dataframe(df, n_bins=3):
    result = pd.DataFrame(index=df.index, columns=df.columns)

    for col in df.columns:
        result[col] = safe_discretize(df[col].values, n_bins)

    return result.values.astype(np.int32)

# 对每个目标基因运行PC发现
def run_ciciv_pc_discovery(data, target_list, method="pc_simple", alpha=0.01):
    pc_results = {}

    print(f"\n  Running {method} for {len(target_list)} target genes...")
    start_time = time.process_time()

    for i, target_idx in enumerate(target_list):
        try:
            if method == "pc_simple":
                pc, ci_num = pc_simple(data, target_idx, alpha, True)
            elif method == "MMPC":
                pc, _, ci_num = MMPC(data, target_idx, alpha, True)
            elif method == "HITON_PC":
                pc, _, ci_num = HITON_PC(data, target_idx, alpha, True)
            else:
                pc = []
            pc_results[target_idx] = pc
        except Exception as e:
            print(f"    Warning: Error for target {target_idx}: {e}")
            pc_results[target_idx] = []

        if (i + 1) % 10 == 0:
            print(f"    Processed {i + 1}/{len(target_list)} genes...")

    end_time = time.process_time()
    return pc_results, end_time - start_time

# 计算调整后的ATE
def calculate_ate_with_adjustment(data, treatment_idx, outcome_idx, adjustment_set):
    treatment = data[:, treatment_idx]
    outcome = data[:, outcome_idx]

    high_mask = treatment == 2
    low_mask = treatment == 0

    if high_mask.sum() == 0 or low_mask.sum() == 0:
        return 0.0

    if len(adjustment_set) > 0 and len(adjustment_set) <= 5:
        ate_sum = 0
        weight_sum = 0

        adj_data = data[:, adjustment_set]
        unique_strata, inverse = np.unique(adj_data, axis=0, return_inverse=True)

        for stratum_idx in range(min(len(unique_strata), 50)):
            stratum_mask = inverse == stratum_idx
            stratum_high = high_mask & stratum_mask
            stratum_low = low_mask & stratum_mask

            if stratum_high.sum() > 0 and stratum_low.sum() > 0:
                p_high = outcome[stratum_high].mean()
                p_low = outcome[stratum_low].mean()
                weight = stratum_mask.sum()
                ate_sum += (p_high - p_low) * weight
                weight_sum += weight

        if weight_sum > 0:
            return ate_sum / weight_sum

    return outcome[high_mask].mean() - outcome[low_mask].mean()

def main():
    TAR_PATH = r"C:\Users\62589\OneDrive\Desktop\brca_metabric.tar.gz"
    OUTPUT_DIR = r"C:\Users\62589\OneDrive\Desktop\Bib-jiaojie\PC-simple\pyCausalFS\CBD\output"
    METHOD = "pc_simple"
    ALPHA = 0.01

    # Step 1: 加载数据
    data_dir = extract_metabric(TAR_PATH)
    expr_df, survival, clinical_df = load_metabric_data(data_dir)

    print(f"  Samples: {len(survival)}")
    print(f"  Total genes: {expr_df.shape[1]}")
    print(f"  Survival events: {survival.sum()} ({survival.mean() * 100:.1f}%)")

    # Step 2: 只使用40个目标基因
    available_genes = [g for g in TARGET_GENES if g in expr_df.columns]
    print(f"  Target genes found: {len(available_genes)}/40")

    # 只取40个目标基因
    expr_subset = expr_df[available_genes].copy()

    print(f"  Subset shape: {expr_subset.shape}")

    # 检查数据
    print(f"  NaN count: {expr_subset.isna().sum().sum()}")
    print(f"  Data range: [{expr_subset.min().min():.2f}, {expr_subset.max().max():.2f}]")

    # Step 3: 离散化
    # 对表达数据离散化
    expr_discrete = discretize_dataframe(expr_subset)

    # 添加SURVIVAL列
    survival_arr = survival.values.astype(np.int32).reshape(-1, 1)
    data_with_survival = np.hstack([expr_discrete, survival_arr])

    gene_cols = list(available_genes)
    survival_idx = len(gene_cols)  # SURVIVAL的列索引

    print(f"  Final data shape: {data_with_survival.shape}")
    print(f"  Discrete value range: [{data_with_survival.min()}, {data_with_survival.max()}]")
    print(f"  SURVIVAL index: {survival_idx}")

    # 验证数据
    assert data_with_survival.min() >= 0, "Data contains negative values!"
    assert data_with_survival.max() <= 2, "Data contains values > 2!"

    # Step 4: PC发现
    target_idx_list = list(range(len(gene_cols)))  # 0到39
    pc_results, runtime = run_ciciv_pc_discovery(
        data_with_survival, target_idx_list, METHOD, ALPHA
    )

    print(f"  Runtime: {runtime:.2f}s")

    # Step 5: 计算ATE
    results = []
    for i, gene in enumerate(available_genes):
        gene_idx = i
        pc_indices = pc_results.get(gene_idx, [])

        # 过滤PC（排除SURVIVAL）
        adjustment_set = [j for j in pc_indices if j != survival_idx and j < len(gene_cols)]

        # 计算ATE
        ate = calculate_ate_with_adjustment(
            data_with_survival, gene_idx, survival_idx, adjustment_set
        )

        # PC基因名称
        pc_gene_names = [gene_cols[j] for j in adjustment_set if j < len(gene_cols)]

        # TCGA对比
        tcga_ate = TCGA_ATE.get(gene, 0)
        tcga_rank = TARGET_GENES.index(gene) + 1
        sign_match = (ate > 0) == (tcga_ate > 0)

        results.append({
            "Gene": gene,
            "Gene_Index": gene_idx,
            "TCGA_Rank": tcga_rank,
            "TCGA_ATE": tcga_ate,
            "METABRIC_ATE": ate,
            "Sign_Match": sign_match,
            "PC_Count": len(pc_indices),
            "PC_Genes": pc_gene_names,
            "Survival_in_PC": survival_idx in pc_indices
        })

    results_df = pd.DataFrame(results)
    results_df["MB_ATE_Rank"] = results_df["METABRIC_ATE"].abs().rank(ascending=False).astype(int)
    results_df = results_df.sort_values("TCGA_Rank")

    # Step 6: 统计指标
    sign_match_rate = results_df["Sign_Match"].mean()
    ate_pearson = np.corrcoef(results_df["TCGA_ATE"], results_df["METABRIC_ATE"])[0, 1]
    ate_spearman, ate_sp_p = stats.spearmanr(results_df["TCGA_ATE"], results_df["METABRIC_ATE"])
    abs_corr, abs_p = stats.spearmanr(results_df["TCGA_ATE"].abs(), results_df["METABRIC_ATE"].abs())
    rank_corr, rank_p = stats.spearmanr(results_df["TCGA_Rank"], results_df["MB_ATE_Rank"])

    top10_tcga = set(results_df.nsmallest(10, "TCGA_Rank")["Gene"])
    top10_mb = set(results_df.nsmallest(10, "MB_ATE_Rank")["Gene"])
    top10_overlap = len(top10_tcga & top10_mb)

    top20_tcga = set(results_df.nsmallest(20, "TCGA_Rank")["Gene"])
    top20_mb = set(results_df.nsmallest(20, "MB_ATE_Rank")["Gene"])
    top20_overlap = len(top20_tcga & top20_mb)

    survival_pc_count = results_df["Survival_in_PC"].sum()

    # 输出
    print(f"\n[A] Dataset")
    print(f"    METABRIC samples: {len(survival)}")
    print(f"    Genes validated: {len(available_genes)}/40")

    print(f"\n[B] ATE Direction Consistency")
    print(
        f"    Sign match: {sign_match_rate * 100:.1f}% ({int(sign_match_rate * len(available_genes))}/{len(available_genes)})")

    print(f"\n[C] ATE Correlation")
    print(f"    Pearson: r = {ate_pearson:.4f}")
    print(f"    Spearman: rho = {ate_spearman:.4f}, p = {ate_sp_p:.4f}")
    print(f"    |ATE| Spearman: rho = {abs_corr:.4f}, p = {abs_p:.4f}")

    print(f"\n[D] Ranking")
    print(f"    Rank correlation: rho = {rank_corr:.4f}, p = {rank_p:.4f}")
    print(f"    Top-10 overlap: {top10_overlap}/10")
    print(f"    Top-20 overlap: {top20_overlap}/20")

    print(f"\n[E] PC Discovery")
    print(f"    Genes with SURVIVAL as PC: {survival_pc_count}")
    print(f"    Avg PC count: {results_df['PC_Count'].mean():.1f}")

    print(f"\n[F] Top 20 Genes")
    print(
        f"{'Rank':<5} {'Gene':<10} {'TCGA_ATE':>10} {'MB_ATE':>10} {'Match':>6} {'MB_Rank':>8} {'#PC':>5} {'Surv':>5}")
    for _, row in results_df.head(20).iterrows():
        match = "✓" if row["Sign_Match"] else "✗"
        surv = "Yes" if row["Survival_in_PC"] else ""
        print(
            f"{row['TCGA_Rank']:<5} {row['Gene']:<10} {row['TCGA_ATE']:>10.4f} {row['METABRIC_ATE']:>10.4f} {match:>6} {row['MB_ATE_Rank']:>8} {row['PC_Count']:>5} {surv:>5}")

    # 保存
    pc_path = os.path.join(OUTPUT_DIR, "metabric_pc.txt")
    with open(pc_path, "w") as f:
        for i, gene in enumerate(available_genes):
            pc = pc_results.get(i, [])
            f.write(f"the pc of {i} ({gene}) is:{pc}\n")
        f.write(f"the running time is: {runtime}\n")

    report_path = os.path.join(OUTPUT_DIR, "metabric_ciciv_report.txt")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("CICIV Framework - METABRIC External Validation\n")

        f.write("1. DATASET\n")
        f.write(f"   Samples: {len(survival)}\n")
        f.write(f"   Genes: {len(available_genes)}/40\n")
        f.write(f"   Method: {METHOD}, Alpha: {ALPHA}\n")
        f.write(f"   Runtime: {runtime:.2f}s\n\n")

        f.write("2. VALIDATION METRICS\n")
        f.write(f"   Sign Consistency: {sign_match_rate * 100:.1f}%\n")
        f.write(f"   Pearson (ATE): r = {ate_pearson:.4f}\n")
        f.write(f"   Spearman (ATE): rho = {ate_spearman:.4f}, p = {ate_sp_p:.4f}\n")
        f.write(f"   Spearman (|ATE|): rho = {abs_corr:.4f}, p = {abs_p:.4f}\n")
        f.write(f"   Rank Correlation: rho = {rank_corr:.4f}, p = {rank_p:.4f}\n")
        f.write(f"   Top-10 Overlap: {top10_overlap}/10\n")
        f.write(f"   Top-20 Overlap: {top20_overlap}/20\n\n")

        f.write("3. PC DISCOVERY\n")
        f.write(f"   Genes with SURVIVAL as PC: {survival_pc_count}\n")
        f.write(f"   Avg PC count: {results_df['PC_Count'].mean():.1f}\n\n")

        f.write("4. DETAILED RESULTS\n")
        f.write(
            f"{'Rank':<5} {'Gene':<10} {'TCGA_ATE':>10} {'MB_ATE':>10} {'Match':>6} {'MB_Rank':>8} {'#PC':>5} {'PC_Genes'}\n")
        for _, row in results_df.iterrows():
            match = "Yes" if row["Sign_Match"] else "No"
            pc_str = ", ".join(row["PC_Genes"][:5]) if row["PC_Genes"] else "-"
            f.write(
                f"{row['TCGA_Rank']:<5} {row['Gene']:<10} {row['TCGA_ATE']:>10.4f} {row['METABRIC_ATE']:>10.4f} {match:>6} {row['MB_ATE_Rank']:>8} {row['PC_Count']:>5} {pc_str}\n")

    csv_path = os.path.join(OUTPUT_DIR, "metabric_ciciv_results.csv")
    results_df.to_csv(csv_path, index=False)

    print(f"Saved: {pc_path}")
    print(f"Saved: {report_path}")
    print(f"Saved: {csv_path}")

if __name__ == "__main__":
    main()