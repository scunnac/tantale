#!/usr/bin/env python3
"""Read PrediTALE's RVD specificities out of its model file (ledger §55).

PrediTALE ships its fitted model as Jstacs XML inside the jar:
  projects/tals/prediction/preditale_quantitative_PBM.xml
Two blocks hold the specificities, both as unnormalised log-scores over
A, C, G, T:
  specsThirteen  21 x 4, keyed on the 13th residue (the RVD's second
                 letter), which is how rare RVDs get a specificity
  separateSpecs   6 x 4, for the RVDs fitted on their own (separateMap
                 names them: HD, NN, NG, HG, NI, NK)
The sibling blocks `a` and `b` (length 905) are optimiser state, not
specificities -- reading them gives values that contradict the biology,
which is how the mistake shows up.

Validated: specsThirteen gives D->C, G->T, I->A, K->G, N->G, and
separateSpecs gives HD->C, NG->T, NI->A, NN->G, NK->G.

Usage:
  python3 dev/preditale-specificities.py <path to PrediTALE.jar> [out.csv]
The jar is GPL-3 (inst/COPYRIGHTS); the table it yields is derived from
that work, so think about licensing before shipping it in tantale.
"""
import csv, math, re, sys, zipfile

MEMBER = "projects/tals/prediction/preditale_quantitative_PBM.xml"
NUM = r'<pos val="\d+">(-?\d+\.?\d*(?:[eE][-+]?\d+)?)</pos>'


def softmax(v):
    m = max(v)
    e = [math.exp(x - m) for x in v]
    t = sum(e)
    return [x / t for x in e]


def main(jar, out="preditale_specificities.csv"):
    with zipfile.ZipFile(jar) as z:
        s = z.read(MEMBER).decode("utf-8", "replace")

    def block(tag):
        i = s.index("<" + tag + ">")
        return s[i:s.index("</" + tag + ">", i)]

    def nums(tag):
        return [float(x) for x in re.findall(NUM, block(tag))]

    thirteen = [x for x in re.findall(r'<pos val="\d+">([^<]+)</pos>',
                                      block("thirteen")) if len(x) == 1]
    symbols = re.findall(
        r'<pos val="\d+">([^<]+)</pos>',
        re.search(r"<symbols>.*?<length>420</length>(.*?)</symbols>", s, re.S).group(1))
    sep_map = [int(x) for x in re.findall(r'<pos val="\d+">(-?\d+)</pos>',
                                          block("separateMap"))]
    separate = {g: r for r, g in zip(symbols, sep_map) if g != -1}

    st, ss = nums("specsThirteen"), nums("separateSpecs")
    rows = [["table", "key", "A", "C", "G", "T", "preferred"]]
    for i, aa in enumerate(thirteen):
        p = softmax(st[4 * i:4 * i + 4])
        rows.append(["by_13th_residue", aa] + [round(x, 4) for x in p] +
                    ["ACGT"[p.index(max(p))]])
    for g in sorted(separate):
        p = softmax(ss[4 * g:4 * g + 4])
        rows.append(["separate_rvd", separate[g]] + [round(x, 4) for x in p] +
                    ["ACGT"[p.index(max(p))]])

    with open(out, "w", newline="") as f:
        csv.writer(f).writerows(rows)
    print(f"{len(rows) - 1} rows -> {out}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(*sys.argv[1:3])
