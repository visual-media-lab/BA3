# coding: utf-8
import argparse
import pandas as pd
import re

def pos_key(mut):
    # 変異表記から最初の数字を抽出してソートキーにする（AA/塩基どちらでもOK）
    m = re.search(r"(\d+)", str(mut))
    return (int(m.group(1)) if m else 10**9, str(mut))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--A", required=True, help="CSV of list A")
    ap.add_argument("--B", required=True, help="CSV of list B")
    ap.add_argument("--out", default="compare_A_B.csv")
    args = ap.parse_args()

    A = pd.read_csv(args.A, dtype=str)
    B = pd.read_csv(args.B, dtype=str)

    a_col = A.columns[0]
    b_col = B.columns[0]

    setA = set(A[a_col].dropna().astype(str))
    setB = set(B[b_col].dropna().astype(str))

    onlyA = sorted(setA - setB, key=pos_key)
    both = sorted(setA & setB, key=pos_key)
    onlyB = sorted(setB - setA, key=pos_key)

    out = pd.DataFrame({
        "Only_A": pd.Series(onlyA),
        "Both": pd.Series(both),
        "Only_B": pd.Series(onlyB)
    })

    out.to_csv(args.out, index=False, encoding="utf-8-sig")
    print("Saved:", args.out)
    print("Only_A:", len(onlyA), "Both:", len(both), "Only_B:", len(onlyB))

if __name__ == "__main__":
    main()