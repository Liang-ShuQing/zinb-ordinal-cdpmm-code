############################################################
# SM tables 5--9: CDPMM (main-paper joint summaries) vs Gaussian
# (HPC gauss-only S=500). Format matches tables 1--4.
############################################################
from __future__ import annotations

import csv
import sys
import re
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

PROJ = Path(__file__).resolve().parents[2]
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else PROJ / "sim_cdpmm_vs_gauss_tables_rest.tex"
BIAS_POOR = Decimal("0.10")
CP_POOR = Decimal("0.90")

MVT = PROJ / "hpc_jobs" / "mvt8_sc12_n200n400_s500_taiyuan_20260915"
S3 = PROJ / "hpc_jobs" / "scen3v2_mvt8_n200n400_s500_wuzhen_20260915"
G12 = PROJ / "hpc_jobs" / "sm_gauss_sc12_mvt8_s500_taiyuan_20260918"
G3 = PROJ / "hpc_jobs" / "sm_gauss_sc3_s500_wuzhen_20260918"

TABLES = [
    (
        "tab:cdpmm_vs_gauss_mvt_s1",
        "multivariate-$t$ random effects, $df=8$",
        1,
        {
            200: (
                MVT / "s1_n200_mvt8" / "joint_param_summary_scen1_mvt_joint_vs_sep_n200_chain10000.csv",
                G12 / "joint_param_summary_gauss_scen1_mvt_gauss_only_n200_chain10000.csv",
            ),
            400: (
                MVT / "s1_n400_mvt8" / "joint_param_summary_scen1_mvt_joint_vs_sep_n400_chain10000.csv",
                G12 / "joint_param_summary_gauss_scen1_mvt_gauss_only_n400_chain10000.csv",
            ),
        },
    ),
    (
        "tab:cdpmm_vs_gauss_mvt_s2",
        "multivariate-$t$ random effects, $df=8$",
        2,
        {
            200: (
                MVT / "s2_n200_mvt8" / "joint_param_summary_scen2_mvt_joint_vs_sep_n200_chain10000.csv",
                G12 / "joint_param_summary_gauss_scen2_mvt_gauss_only_n200_chain10000.csv",
            ),
            400: (
                MVT / "s2_n400_mvt8" / "joint_param_summary_scen2_mvt_joint_vs_sep_n400_chain10000.csv",
                G12 / "joint_param_summary_gauss_scen2_mvt_gauss_only_n400_chain10000.csv",
            ),
        },
    ),
    (
        "tab:cdpmm_vs_gauss_normal_s3",
        "multivariate normal random effects",
        3,
        {
            200: (
                S3 / "s3_n200_normal" / "joint_param_summary_scen3_normal_joint_vs_sep_n200_chain10000.csv",
                G3 / "joint_param_summary_gauss_scen3_normal_gauss_only_n200_chain10000.csv",
            ),
            400: (
                S3 / "s3_n400_normal" / "joint_param_summary_scen3_normal_joint_vs_sep_n400_chain10000.csv",
                G3 / "joint_param_summary_gauss_scen3_normal_gauss_only_n400_chain10000.csv",
            ),
        },
    ),
    (
        "tab:cdpmm_vs_gauss_mixture_s3",
        "multivariate mixture random effects",
        3,
        {
            200: (
                S3 / "s3_n200_mixture" / "joint_param_summary_scen3_mixture_joint_vs_sep_n200_chain10000.csv",
                G3 / "joint_param_summary_gauss_scen3_mixture_gauss_only_n200_chain10000.csv",
            ),
            400: (
                S3 / "s3_n400_mixture" / "joint_param_summary_scen3_mixture_joint_vs_sep_n400_chain10000.csv",
                G3 / "joint_param_summary_gauss_scen3_mixture_gauss_only_n400_chain10000.csv",
            ),
        },
    ),
    (
        "tab:cdpmm_vs_gauss_mvt_s3",
        "multivariate-$t$ random effects, $df=8$",
        3,
        {
            200: (
                S3 / "s3_n200_mvt8" / "joint_param_summary_scen3_mvt_joint_vs_sep_n200_chain10000.csv",
                G3 / "joint_param_summary_gauss_scen3_mvt_gauss_only_n200_chain10000.csv",
            ),
            400: (
                S3 / "s3_n400_mvt8" / "joint_param_summary_scen3_mvt_joint_vs_sep_n400_chain10000.csv",
                G3 / "joint_param_summary_gauss_scen3_mvt_gauss_only_n400_chain10000.csv",
            ),
        },
    ),
]


def r3(x: float) -> Decimal:
    return Decimal(str(x)).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP)


def fmt_true(x: float) -> str:
    d = Decimal(str(x))
    one = d.quantize(Decimal("0.1"), rounding=ROUND_HALF_UP)
    if abs(d - one) < Decimal("0.0000001"):
        s = f"{one:.1f}"
    else:
        s = f"{r3(x):.3f}"
    if s.startswith("-"):
        return f"$-{s[1:]}$"
    return s


def fmt_num(x: float, *, bias: bool = False, cp: bool = False) -> str:
    d = r3(x)
    s = f"{d:.3f}"
    if s.startswith("-"):
        tex = f"$-{s[1:]}$"
    else:
        tex = s
    bold = (bias and abs(d) >= BIAS_POOR) or (cp and d < CP_POOR)
    if bold:
        return f"\\textbf{{{tex}}}"
    return tex


def tex_param(name: str) -> str:
    m = re.match(r"^(alpha|beta|gamma)(\d+)$", name)
    if m:
        sym = {"alpha": r"\alpha", "beta": r"\beta", "gamma": r"\gamma"}[m.group(1)]
        return f"${sym}_{{{m.group(2)}}}$"
    if name == "r":
        return "$r$"
    m = re.match(r"^Sigma([123])\1$", name)
    if m:
        return f"$\\sigma_{m.group(1)}^2$"
    m = re.match(r"^rho([123])$", name)
    if m:
        pair = {"1": "12", "2": "13", "3": "23"}[m.group(1)]
        return f"$\\rho_{{{pair}}}$"
    raise ValueError(name)


def load(path: Path) -> dict[str, dict]:
    if not path.is_file():
        raise FileNotFoundError(path)
    rows = {}
    with path.open(encoding="utf-8-sig", newline="") as f:
        for rec in csv.DictReader(f):
            rows[rec["Parameter"]] = rec
    return rows


def keep(name: str) -> bool:
    return not name.startswith("delta")


def make_table(label: str, re_lab: str, scen: int, blocks: dict) -> list[str]:
    cap = (
        f"Bias, RMSE, and coverage under a CDPMM versus Gaussian joint prior "
        f"({re_lab}; Scenario~{scen}, $N\\in\\{{200,400\\}}$). "
        f"Bold Bias entries mark $|\\mathrm{{Bias}}|\\ge 0.10$; "
        f"bold CP entries mark $\\mathrm{{CP}}<0.90$ (RMSE is not bolded)."
    )
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        rf"\caption{{{cap}}}",
        rf"\label{{{label}}}",
        r"\begin{adjustbox}{max width=\textwidth,center}",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{2.8pt}",
        r"\begin{tabular}{lccccccc}",
        r"\toprule",
        r" & & \multicolumn{3}{c}{\textbf{CDPMM Joint}} & \multicolumn{3}{c}{\textbf{Gaussian Joint}} \\",
        r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}",
        r"\textbf{Parameter} & \textbf{True Value} & Bias & RMSE & CP & Bias & RMSE & CP \\",
        r"\midrule",
    ]
    for n in (200, 400):
        cdp, gau = blocks[n]
        cd = load(cdp)
        ga = load(gau)
        names = [k for k in cd if keep(k)]
        missing = [k for k in names if k not in ga]
        if missing:
            raise SystemExit(f"missing in gauss {gau.name}: {missing}")
        lines.append(rf"\multicolumn{{8}}{{@{{}}l}}{{\textbf{{$N={n}$}}}} \\")
        lines.append(r"\midrule")
        for k in names:
            tv = float(cd[k]["True_Value"])
            row = " & ".join(
                [
                    tex_param(k),
                    fmt_true(tv),
                    fmt_num(float(cd[k]["Bias"]), bias=True),
                    fmt_num(float(cd[k]["RMSE"])),
                    fmt_num(float(cd[k]["CP"]), cp=True),
                    fmt_num(float(ga[k]["Bias"]), bias=True),
                    fmt_num(float(ga[k]["RMSE"])),
                    fmt_num(float(ga[k]["CP"]), cp=True),
                ]
            )
            lines.append(row + r" \\")
        if n == 200:
            lines.append(r"\midrule")
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{adjustbox}", r"\end{table}", ""])
    return lines


def main():
    out = [
        "% Tables 5--9: CDPMM from main-paper joint summaries; Gaussian from HPC gauss-only S=500.",
        "% Same columns and bolding as sim_cdpmm_vs_gauss_tables.tex (tables 1--4).",
        "",
    ]
    for label, re_lab, scen, blocks in TABLES:
        out.extend(make_table(label, re_lab, scen, blocks))
        print("ok", label)
    OUT.write_text("\n".join(out) + "\n", encoding="utf-8")
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
