# coding: utf-8
import argparse
import pandas as pd
import re
from collections import Counter

def parse_list_field(x):
    if pd.isna(x) or str(x).strip() == "":
        return []
    return [s.strip() for s in str(x).split(",") if s.strip()]

def extract_spike_aa(aa_subs):
    # Nextclade aaSubstitutions: "S:A67V,ORF1a:T3255I,..."
    muts = []
    for s in parse_list_field(aa_subs):
        if s.startswith("S:"):
            muts.append(s[2:])
    return muts

def get_accession_from_seqname(seqname):
    m = re.search(r"\|(EPI_ISL_\d+)\|", str(seqname))
    return m.group(1) if m else ""

def year_filter(df, start_date, end_date):
    df["Collection date"] = df["Collection date"].astype(str)
    # YYYY-MM-DD の先頭10文字だけ比較してもOK
    df = df[(df["Collection date"] >= start_date) & (df["Collection date"] <= end_date)]
    return df

def common_muts(list_series, threshold=0.6):
    """
    list_series: 変異リストのSeries（各行がlist）
    threshold: 出現率の下限（0.5なら50%以上）
    """
    n = len(list_series)
    c = Counter()
    for muts in list_series:
        c.update(set(muts))  # 同一サンプル内の重複は1回だけ

    out = []
    for m, count in c.items():
        freq = count / n
        if freq >= threshold:
            out.append((m, count, freq))

    out.sort(key=lambda x: (-x[2], x[0]))
    return out, n

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nextclade", required=True, help="nextclade.tsv (or csv/tsv). Must contain seqName, aaSubstitutions, substitutions")
    ap.add_argument("--meta", required=True, help="GISAID metadata.tsv containing Accession ID, Collection date")
    ap.add_argument("--start", required=True, help="start date YYYY-MM-DD")
    ap.add_argument("--end", required=True, help="end date YYYY-MM-DD")
    ap.add_argument("--threshold", type=float, default=0.6, help="frequency threshold (default 0.5)")
    args = ap.parse_args()

    # Nextclade読み込み（tsv/csv自動判別）
    if args.nextclade.lower().endswith(".tsv"):
        nc = pd.read_csv(args.nextclade, sep="\t", dtype=str)
    else:
        nc = pd.read_csv(args.nextclade, dtype=str)

    meta = pd.read_csv(args.meta, sep="\t", dtype=str)

    # Accession ID を合わせる
    if "Accession ID" not in nc.columns:
        nc["Accession ID"] = nc["seqName"].apply(get_accession_from_seqname)

    df = nc.merge(meta, on="Accession ID", how="left")
    df = df.dropna(subset=["Collection date"])
    df = year_filter(df, args.start, args.end)

    print(f"Records in range {args.start} - {args.end}:", len(df))

    # 変異をlistにする
    df["NuclSubs"] = df["substitutions"].apply(parse_list_field)
    df["SpikeAA"] = df["aaSubstitutions"].apply(extract_spike_aa)

    # 共通変異を抽出
    common_nucl, n1 = common_muts(df["NuclSubs"], threshold=args.threshold)
    common_spike, n2 = common_muts(df["SpikeAA"], threshold=args.threshold)

    out_nucl = pd.DataFrame(common_nucl, columns=["Nucl_mut", "count", "freq"])
    out_spike = pd.DataFrame(common_spike, columns=["Spike_AA_mut", "count", "freq"])

    # 保存
    tag = f"{args.start}_to_{args.end}".replace("-", "")
    out_nucl.to_csv(f"common_nucl_{tag}.csv", index=False, encoding="utf-8-sig")
    out_spike.to_csv(f"common_spikeaa_{tag}.csv", index=False, encoding="utf-8-sig")

    print("Saved:")
    print(f" common_nucl_{tag}.csv ({len(out_nucl)})")
    print(f" common_spikeaa_{tag}.csv ({len(out_spike)})")

if __name__ == "__main__":
    main()