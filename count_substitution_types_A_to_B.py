# coding: utf-8
import argparse
import pandas as pd
import re

TARGETS = [
    "C>U", "G>A", "A>G", "U>C", "G>U", "C>A",
    "G>C", "C>G", "A>C", "U>G", "A>U", "U>A"
]

def normalize_nt(nt: str) -> str:
    nt = nt.upper()
    return "U" if nt == "T" else nt

def parse_point_substitution(mut: str):
    """
    Parse simple nucleotide substitution like C193T or A23403G.
    Returns tuple (ref_base, alt_base) as RNA bases (A,C,G,U),
    or None if not a simple point substitution.
    """
    if mut is None:
        return None

    s = str(mut).strip().upper()

    # common cleanup (e.g., trailing '+', spaces)
    s = s.strip().rstrip("+")

    # Must match: <base><pos><base>
    m = re.match(r"^([ACGT])(\d+)([ACGT])$", s)
    if not m:
        return None

    ref = normalize_nt(m.group(1))
    alt = normalize_nt(m.group(3))
    return ref, alt

def to_label(ref_alt):
    ref, alt = ref_alt
    return f"{ref}>{alt}"

def reverse_ref_alt(ref_alt):
    ref, alt = ref_alt
    return (alt, ref)

def count_list(muts, reverse=False):
    """
    Count substitution types in a list of mutation strings.
    reverse=True means count reverse direction (alt->ref).
    Returns dict counts + other count.
    """
    counts = {k: 0 for k in TARGETS}
    other = 0

    for mut in muts:
        ra = parse_point_substitution(mut)
        if ra is None:
            other += 1
            continue

        if reverse:
            ra = reverse_ref_alt(ra)

        label = to_label(ra)
        if label in counts:
            counts[label] += 1
        else:
            other += 1

    return counts, other

def main():
    ap = argparse.ArgumentParser(
        description="Compute A->B mutation spectrum from compare_mut_lists output. "
                    "Only_B counted forward; Only_A counted in reverse; Both ignored."
    )
    ap.add_argument("--compare", required=True, help="compare_A_B.csv (columns: Only_A, Both, Only_B)")
    ap.add_argument("--out", default="mutation_spectrum_A_to_B.csv")
    args = ap.parse_args()

    df = pd.read_csv(args.compare, dtype=str)

    # tolerate slightly different column names
    colA = "Only_A" if "Only_A" in df.columns else df.columns[0]
    colB = "Only_B" if "Only_B" in df.columns else df.columns[-1]

    only_a = df[colA].dropna().astype(str).tolist()  # present in A but not in B -> reverse direction
    only_b = df[colB].dropna().astype(str).tolist()  # present in B but not in A -> forward direction

    forward_counts, forward_other = count_list(only_b, reverse=False)
    reverse_counts, reverse_other = count_list(only_a, reverse=True)

    # Sum them -> A->B spectrum
    total = {k: forward_counts[k] + reverse_counts[k] for k in TARGETS}
    total_other = forward_other + reverse_other

    out_rows = []
    for k in TARGETS:
        out_rows.append({
            "Substitution": k,
            "Count_total(A_to_B)": total[k],
            "Count_forward(Only_B)": forward_counts[k],
            "Count_reverse(Only_A)": reverse_counts[k],
        })

    out_rows.append({
        "Substitution": "OTHER/Non-simple",
        "Count_total(A_to_B)": total_other,
        "Count_forward(Only_B)": forward_other,
        "Count_reverse(Only_A)": reverse_other,
    })

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(args.out, index=False, encoding="utf-8-sig")

    print("Saved:", args.out)
    print(out_df.to_string(index=False))

if __name__ == "__main__":
    main()