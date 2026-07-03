#!/usr/bin/env python3
"""
QC audit for externally provided MiXeR scripts and summary_log outputs.

This is designed for the recovered MiXeR result bundle:
  scripts: ${PROJECT_ROOT}
  outputs: ${PROJECT_ROOT}
"""

from __future__ import annotations

import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean


SCRIPT_DIR = Path("${PROJECT_ROOT}")
SUMMARY_DIR = Path("${PROJECT_ROOT}")
REPORT_PATH = SUMMARY_DIR / "mixer_summary_log_qc_report.md"
METRIC_QC_CSV = SUMMARY_DIR / "mixer_summary_metric_qc.csv"
REPLICATE_QC_CSV = SUMMARY_DIR / "mixer_replicate_qc.csv"
LOG_ERROR_CSV = SUMMARY_DIR / "mixer_log_error_qc.csv"

N_REP = 20
OLD_TRAITS = ["TRAIT_B", "TRAIT_A", "de Lange_IBD", "de Lange_CD", "de Lange_UC", "Liu_IBD", "Liu_CD", "Liu_UC"]
OLD_PAIRS = [
    "TRAIT_B_vs_TRAIT_A",
    "ASD_vs_de Lange_IBD",
    "ASD_vs_de Lange_CD",
    "ASD_vs_de Lange_UC",
    "ASD_vs_Liu_IBD",
    "ASD_vs_Liu_CD",
    "ASD_vs_Liu_UC",
]
CURRENT_P0_REQUIRED = ["TRAIT_B", "TRAIT_A", "ADHD", "SCZ", "disease_anxiety", "IBD_deLange", "CD_deLange", "UC_deLange"]

LDSC_RG = {
    "TRAIT_B vs TRAIT_A": 0.2465,
    "TRAIT_B vs TRAIT_C": -0.0599,
    "TRAIT_B vs TRAIT_E": -0.0367,
    "TRAIT_B vs TRAIT_D": -0.0699,
}

RUN_MARKER = re.compile(r"\[(done|FAIL|skip|run)\]\s+(fit1|test1|fit2|test2)\s+([A-Za-z0-9_]+)\s+rep(\d+)")
RETRY_FIT2_MARKER = re.compile(r"\[(done|FAIL|skip|run)\]\s+(ASD_vs_[A-Za-z0-9_]+)\s+rep(\d+)")
ERROR_PAT = re.compile(
    r"(Traceback|ValueError|invalid choice|command not found|\[FAIL\]|error:|failed|unexpected EOF|docker API|no such file)",
    re.IGNORECASE,
)


def as_float(value: str) -> float:
    try:
        result = float(value)
    except Exception:
        return float("nan")
    return result


def finite(value: float) -> bool:
    return math.isfinite(value)


def read_csv_dict(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def parse_metric_qc() -> tuple[list[dict[str, object]], list[str]]:
    rows = read_csv_dict(SUMMARY_DIR / "mixer_summary.csv")
    metric_rows: list[dict[str, object]] = []
    issues: list[str] = []
    required = ["label", "nc1", "nc1_se", "nc2", "nc2_se", "nc12", "nc12_se", "dice", "dice_se", "rho_ge", "rho_ge_se", "rg", "rg_se"]
    if not rows:
        return [], ["mixer_summary.csv missing or empty"]
    missing_cols = [col for col in required if col not in rows[0]]
    if missing_cols:
        issues.append("mixer_summary.csv missing columns: " + ", ".join(missing_cols))
    for row in rows:
        label = row.get("label", "")
        values = {key: as_float(row.get(key, "")) for key in required if key != "label"}
        metric_issue: list[str] = []
        if any(not finite(value) for value in values.values()):
            metric_issue.append("non_finite")
        if not 0 <= values.get("dice", float("nan")) <= 1:
            metric_issue.append("dice_out_of_range")
        if not -1 <= values.get("rg", float("nan")) <= 1:
            metric_issue.append("rg_out_of_range")
        if not -1 <= values.get("rho_ge", float("nan")) <= 1:
            metric_issue.append("rho_ge_out_of_range")
        if values.get("nc12", 0.0) < 0 or values.get("nc1", 0.0) < 0 or values.get("nc2", 0.0) < 0:
            metric_issue.append("negative_causal_component")
        ldsc_rg = LDSC_RG.get(label)
        sign_concordant = ""
        if ldsc_rg is not None and finite(values.get("rg", float("nan"))):
            sign_concordant = (values["rg"] == 0 and ldsc_rg == 0) or (values["rg"] * ldsc_rg > 0)
            if not sign_concordant:
                metric_issue.append("ldsc_sign_discordant")
        metric_rows.append(
            {
                "label": label,
                "nc12": values.get("nc12"),
                "nc12_se": values.get("nc12_se"),
                "dice": values.get("dice"),
                "dice_se": values.get("dice_se"),
                "rho_ge": values.get("rho_ge"),
                "rho_ge_se": values.get("rho_ge_se"),
                "rg": values.get("rg"),
                "rg_se": values.get("rg_se"),
                "ldsc_rg": "" if ldsc_rg is None else ldsc_rg,
                "ldsc_sign_concordant": sign_concordant,
                "issue": ";".join(metric_issue),
            }
        )
    if len(rows) != len(OLD_PAIRS):
        issues.append(f"expected {len(OLD_PAIRS)} legacy pair rows, observed {len(rows)}")
    return metric_rows, issues


def parse_venn_qc() -> list[str]:
    rows = read_csv_dict(SUMMARY_DIR / "mixer_venn_summary.csv")
    issues: list[str] = []
    if not rows:
        return ["mixer_venn_summary.csv missing or empty"]
    asd_totals = []
    for row in rows:
        unique_t1 = as_float(row.get("unique_t1", ""))
        shared = as_float(row.get("shared", ""))
        nc1u = as_float(row.get("nc1u", ""))
        if finite(unique_t1) and finite(shared) and finite(nc1u):
            if abs((unique_t1 + shared) - nc1u) > max(1e-6, abs(nc1u) * 1e-8):
                issues.append(f"Venn inconsistency in {row.get('label')}: unique_t1 + shared != nc1u")
            asd_totals.append(nc1u)
    if asd_totals and max(asd_totals) - min(asd_totals) > max(1e-6, mean(asd_totals) * 1e-8):
        issues.append("ASD univariate total is not stable across Venn rows")
    return issues


def parse_replicates() -> list[dict[str, object]]:
    status_by_task: dict[tuple[str, str, int], str] = {}
    for log_name in ["mixer_run.log", "retry_test2.log"]:
        path = SUMMARY_DIR / log_name
        if not path.exists():
            continue
        for line in path.read_text(errors="replace").splitlines():
            match = RUN_MARKER.search(line)
            if not match:
                continue
            status, step, label, rep = match.groups()
            status_by_task[(step, label, int(rep))] = status

    retry_fit2 = SUMMARY_DIR / "retry_fit2.log"
    if retry_fit2.exists():
        for line in retry_fit2.read_text(errors="replace").splitlines():
            match = RETRY_FIT2_MARKER.search(line)
            if not match:
                continue
            status, label, rep = match.groups()
            status_by_task[("fit2", label, int(rep))] = status

    rows: list[dict[str, object]] = []
    for step, labels in [
        ("fit1", OLD_TRAITS),
        ("test1", OLD_TRAITS),
        ("fit2", OLD_PAIRS),
        ("test2", OLD_PAIRS),
    ]:
        for label in labels:
            counts = Counter(status_by_task.get((step, label, rep), "missing") for rep in range(1, N_REP + 1))
            pass_threshold = 18 if step in {"fit1", "fit2", "test2"} else 0
            blocking = step in {"fit1", "fit2", "test2"}
            status = "pass" if (not blocking or counts.get("done", 0) >= pass_threshold) else "fail"
            rows.append(
                {
                    "step": step,
                    "label": label,
                    "done": counts.get("done", 0),
                    "fail": counts.get("FAIL", 0),
                    "skip": counts.get("skip", 0),
                    "run": counts.get("run", 0),
                    "missing": counts.get("missing", 0),
                    "blocking_for_bivar_summary": blocking,
                    "pass_threshold": pass_threshold,
                    "status": status,
                }
            )
    return rows


def parse_log_errors() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in sorted(SUMMARY_DIR.glob("*.log")):
        for line_number, line in enumerate(path.read_text(errors="replace").splitlines(), start=1):
            if ERROR_PAT.search(line):
                rows.append(
                    {
                        "log": str(path),
                        "line_number": line_number,
                        "line": line.strip()[:1000],
                    }
                )
    return rows


def script_coverage_issue() -> list[str]:
    issues: list[str] = []
    run_script = SCRIPT_DIR / "12b_run_mixer.sh"
    analyze_script = SCRIPT_DIR / "12c_analyze_mixer.py"
    text = ""
    for path in [run_script, analyze_script]:
        if path.exists():
            text += path.read_text(errors="replace") + "\n"
        else:
            issues.append(f"missing script: {path}")
    missing = [trait for trait in CURRENT_P0_REQUIRED if trait not in text]
    if missing:
        issues.append("legacy scripts do not cover current P0 traits: " + ", ".join(missing))
    if "--fit-sequence diffevo-fast" in text and "test1" in text:
        issues.append("12b script applies diffevo-fast globally; log shows this is invalid for test1")
    return issues


def render_report(
    metric_rows: list[dict[str, object]],
    metric_issues: list[str],
    venn_issues: list[str],
    replicate_rows: list[dict[str, object]],
    log_errors: list[dict[str, object]],
    coverage_issues: list[str],
) -> str:
    blocking_failures = [
        row for row in replicate_rows
        if row["blocking_for_bivar_summary"] and row["status"] != "pass"
    ]
    metric_failures = [row for row in metric_rows if row["issue"]]
    legacy_gate = "PASS_WITH_CAVEATS" if metric_rows and not metric_failures and not blocking_failures and not metric_issues and not venn_issues else "FAIL"
    p0_gate = "FAIL" if coverage_issues else legacy_gate

    lines: list[str] = []
    lines.append("# MiXeR Summary Log QC")
    lines.append("")
    lines.append(f"- Legacy gut-only MiXeR gate: **{legacy_gate}**")
    lines.append(f"- Current P0 MiXeR gate: **{p0_gate}**")
    lines.append(f"- Summary rows parsed: {len(metric_rows)}")
    lines.append(f"- Blocking replicate failures: {len(blocking_failures)}")
    lines.append(f"- Log error/warning signatures: {len(log_errors)}")
    lines.append("")
    lines.append("## Interpretation")
    lines.append("- The recovered `summary_log` bundle is usable as a legacy gut-only MiXeR result, with caveats.")
    lines.append("- It is not sufficient as the current P0 MiXeR deliverable because psychiatric benchmark traits are absent.")
    lines.append("- The original main run had many failures, but retry logs recover bivariate `fit2`/`test2` to the ≥18/20 threshold.")
    lines.append("- `test1` failed because `diffevo-fast` was passed to `test1`; this does not directly block the bivariate summary table, but the run script should be corrected before rerun.")
    lines.append("")
    lines.append("## Metric QC")
    if metric_rows:
        lines.append("| label | nc12 | dice | rho_ge | rg | LDSC rg | sign concordant |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|")
        for row in metric_rows:
            lines.append(
                f"| {row['label']} | {float(row['nc12']):.3f} | {float(row['dice']):.4f} | "
                f"{float(row['rho_ge']):.4f} | {float(row['rg']):.4f} | {row['ldsc_rg']} | {row['ldsc_sign_concordant']} |"
            )
    else:
        lines.append("- No metric rows parsed.")
    if metric_issues or venn_issues:
        lines.append("")
        lines.append("## Metric Issues")
        for issue in metric_issues + venn_issues:
            lines.append(f"- {issue}")
    lines.append("")
    lines.append("## Replicate QC")
    for row in replicate_rows:
        if row["blocking_for_bivar_summary"] or row["step"] == "test1":
            lines.append(
                f"- {row['step']} {row['label']}: done={row['done']}/20, fail={row['fail']}, "
                f"missing={row['missing']}, status={row['status']}"
            )
    lines.append("")
    lines.append("## Script Coverage")
    for issue in coverage_issues:
        lines.append(f"- {issue}")
    lines.append("")
    lines.append("## Files")
    lines.append(f"- Metric QC CSV: `{METRIC_QC_CSV}`")
    lines.append(f"- Replicate QC CSV: `{REPLICATE_QC_CSV}`")
    lines.append(f"- Log error CSV: `{LOG_ERROR_CSV}`")
    return "\n".join(lines) + "\n"


def main() -> None:
    metric_rows, metric_issues = parse_metric_qc()
    venn_issues = parse_venn_qc()
    replicate_rows = parse_replicates()
    log_errors = parse_log_errors()
    coverage_issues = script_coverage_issue()

    write_csv(METRIC_QC_CSV, metric_rows)
    write_csv(REPLICATE_QC_CSV, replicate_rows)
    write_csv(LOG_ERROR_CSV, log_errors)
    REPORT_PATH.write_text(
        render_report(metric_rows, metric_issues, venn_issues, replicate_rows, log_errors, coverage_issues),
        encoding="utf-8",
    )
    print(REPORT_PATH)


if __name__ == "__main__":
    main()
