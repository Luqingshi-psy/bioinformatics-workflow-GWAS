#!/usr/bin/env python3
"""
MiXeR QC for the trait pair (A x B) project.

This script audits existing MiXeR inputs, logs, and expected outputs before
using MiXeR results downstream. It is intentionally read-mostly: the only
files it writes are QC tables and a Markdown report.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any


DESKTOP_MIXER_ROOT = Path("${PROJECT_ROOT}")
T7_MIXER_ROOT = Path("${PROJECT_ROOT}")
NEW_RESULTS_ROOT = Path("${PROJECT_ROOT}")
OUT_DIR = NEW_RESULTS_ROOT / "mixer_qc"
N_REP = 20

OLD_UNIV_TRAITS = [
    "TRAIT_B",
    "TRAIT_A",
    "de Lange_IBD",
    "de Lange_CD",
    "de Lange_UC",
    "Liu_IBD",
    "Liu_CD",
    "Liu_UC",
]
OLD_BIVAR_LABELS = [
    "TRAIT_B_vs_TRAIT_A",
    "ASD_vs_de Lange_IBD",
    "ASD_vs_de Lange_CD",
    "ASD_vs_de Lange_UC",
    "ASD_vs_Liu_IBD",
    "ASD_vs_Liu_CD",
    "ASD_vs_Liu_UC",
]

P0_CURRENT_TRAITS = [
    "TRAIT_B",
    "TRAIT_A",
    "ADHD",
    "SCZ",
    "disease_anxiety",
    "IBD_deLange",
    "CD_deLange",
    "UC_deLange",
]
P0_CURRENT_BIVAR_LABELS = [
    "TRAIT_B_vs_TRAIT_A",
    "ASD_vs_ADHD",
    "ASD_vs_SCZ",
    "ASD_vs_anxiety",
    "TRAIT_B_vs_TRAIT_C",
    "TRAIT_B_vs_TRAIT_E",
    "TRAIT_B_vs_TRAIT_D",
]

ERROR_PAT = re.compile(
    r"(error|failed|traceback|exception|unexpected eof|no such file|cannot connect|docker api|docker.sock)",
    re.IGNORECASE,
)


@dataclass
class OutputExpectation:
    root: str
    label: str
    analysis_type: str
    phase: str
    expected_reps: int
    present_reps: int
    combined_present: bool
    nonempty_reps: int
    status: str


def safe_read_json(path: Path) -> dict[str, Any] | None:
    try:
        with path.open("rt") as handle:
            data = json.load(handle)
        if isinstance(data, dict):
            return data
    except Exception:
        return None
    return None


def get_nested(data: dict[str, Any], keys: list[str]) -> Any:
    current: Any = data
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return None
        current = current[key]
    return current


def finite_or_none(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def ci_mean(data: dict[str, Any], metric: str) -> float | None:
    return finite_or_none(get_nested(data, ["ci", metric, "mean"]))


def audit_expected_outputs(root: Path) -> list[OutputExpectation]:
    rows: list[OutputExpectation] = []
    results_dir = root / "results"
    for label in OLD_UNIV_TRAITS:
        for phase in ["fit", "test"]:
            reps = [results_dir / f"{label}.{phase}.rep{rep}.json" for rep in range(1, N_REP + 1)]
            present = [path for path in reps if path.exists()]
            nonempty = [path for path in present if path.stat().st_size > 0]
            combined = results_dir / f"{label}.{phase}.json"
            status = "pass" if len(nonempty) >= 18 and combined.exists() and combined.stat().st_size > 0 else "fail"
            rows.append(
                OutputExpectation(
                    str(root),
                    label,
                    "univariate",
                    phase,
                    N_REP,
                    len(present),
                    combined.exists() and combined.stat().st_size > 0,
                    len(nonempty),
                    status,
                )
            )
    for label in OLD_BIVAR_LABELS:
        for phase in ["fit", "test"]:
            reps = [results_dir / f"{label}.{phase}.rep{rep}.json" for rep in range(1, N_REP + 1)]
            present = [path for path in reps if path.exists()]
            nonempty = [path for path in present if path.stat().st_size > 0]
            combined = results_dir / f"{label}.{phase}.json"
            status = "pass" if len(nonempty) >= 18 and combined.exists() and combined.stat().st_size > 0 else "fail"
            rows.append(
                OutputExpectation(
                    str(root),
                    label,
                    "bivariate",
                    phase,
                    N_REP,
                    len(present),
                    combined.exists() and combined.stat().st_size > 0,
                    len(nonempty),
                    status,
                )
            )
    return rows


def inventory_sumstats(root: Path) -> list[dict[str, Any]]:
    sumstats_dir = root / "sumstats"
    rows: list[dict[str, Any]] = []
    if not sumstats_dir.exists():
        return rows
    for path in sorted(sumstats_dir.glob("*.sumstats.gz")):
        header = ""
        n_rows = 0
        sample_bad = ""
        try:
            with gzip.open(path, "rt") as handle:
                header = handle.readline().strip()
                for n_rows, line in enumerate(handle, start=1):
                    if n_rows <= 5:
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) != len(header.split("\t")):
                            sample_bad = f"row_{n_rows}_field_count={len(fields)}"
                            break
        except Exception as exc:
            sample_bad = f"read_error={type(exc).__name__}: {exc}"
        columns = header.split("\t") if header else []
        rows.append(
            {
                "root": str(root),
                "file": str(path),
                "trait": path.name.replace(".sumstats.gz", ""),
                "size_mb": round(path.stat().st_size / 1024 / 1024, 3),
                "rows_counted_until_error": n_rows,
                "header": header,
                "has_required_columns": all(col in columns for col in ["SNP", "CHR", "BP", "A1", "A2", "N", "Z"]),
                "sample_issue": sample_bad,
            }
        )
    return rows


def parse_logs(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(root.glob("*.log")):
        try:
            lines = path.read_text(errors="replace").splitlines()
        except Exception as exc:
            rows.append(
                {
                    "root": str(root),
                    "log": str(path),
                    "line_number": 0,
                    "category": "read_error",
                    "line": f"{type(exc).__name__}: {exc}",
                }
            )
            continue
        for idx, line in enumerate(lines, start=1):
            if ERROR_PAT.search(line):
                category = "docker" if "docker" in line.lower() or "container" in line.lower() else "error"
                rows.append(
                    {
                        "root": str(root),
                        "log": str(path),
                        "line_number": idx,
                        "category": category,
                        "line": line.strip()[:1000],
                    }
                )
    return rows


def parse_combined_jsons(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    results_dir = root / "results"
    if not results_dir.exists():
        return rows
    for path in sorted(results_dir.glob("*.json")):
        if ".rep" in path.name:
            continue
        data = safe_read_json(path)
        if data is None:
            rows.append({"root": str(root), "file": str(path), "status": "unreadable"})
            continue
        rows.append(
            {
                "root": str(root),
                "file": str(path),
                "status": "parsed",
                "n1": ci_mean(data, "n1"),
                "n2": ci_mean(data, "n2"),
                "n12": ci_mean(data, "n12"),
                "dice": ci_mean(data, "dice"),
                "rho": ci_mean(data, "rho"),
                "loglike": finite_or_none(data.get("loglike")),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]] | list[OutputExpectation]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    dict_rows: list[dict[str, Any]] = [
        row.__dict__ if hasattr(row, "__dict__") else row for row in rows
    ]
    if not dict_rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(dict_rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(dict_rows)


def summarize_missing_current_traits(sumstats_rows: list[dict[str, Any]]) -> list[str]:
    present = {str(row["trait"]) for row in sumstats_rows}
    return [trait for trait in P0_CURRENT_TRAITS if trait not in present]


def render_report(
    expectation_rows: list[OutputExpectation],
    sumstats_rows: list[dict[str, Any]],
    log_rows: list[dict[str, Any]],
    json_rows: list[dict[str, Any]],
) -> str:
    status_counts = Counter(row.status for row in expectation_rows)
    error_counts = Counter(row["category"] for row in log_rows)
    roots = [DESKTOP_MIXER_ROOT, T7_MIXER_ROOT]
    existing_roots = [str(root) for root in roots if root.exists()]
    missing_current = summarize_missing_current_traits(sumstats_rows)
    combined_parsed = [row for row in json_rows if row.get("status") == "parsed"]
    any_expect_pass = any(row.status == "pass" for row in expectation_rows)
    old_results_dirs = [root / "results" for root in roots]
    old_results_present = [str(path) for path in old_results_dirs if path.exists()]

    gate = "PASS" if any_expect_pass and combined_parsed and not log_rows else "FAIL"

    lines: list[str] = []
    lines.append("# MiXeR QC Report")
    lines.append("")
    lines.append(f"- QC gate: **{gate}**")
    lines.append(f"- Checked roots: {', '.join(existing_roots) if existing_roots else 'none'}")
    lines.append(f"- Results directories found: {', '.join(old_results_present) if old_results_present else 'none'}")
    lines.append(f"- Expected old-design output rows passing replicate gate: {status_counts.get('pass', 0)} / {len(expectation_rows)}")
    lines.append(f"- Combined JSON files parsed: {len(combined_parsed)}")
    lines.append(f"- Log error lines detected: {len(log_rows)}")
    lines.append("")
    lines.append("## Key Findings")
    if not old_results_present:
        lines.append("- No MiXeR `results` directory was found under the checked roots, so fit/test JSON outputs are unavailable.")
    if not combined_parsed:
        lines.append("- No combined MiXeR JSON was parsed; overlap metrics such as n1/n2/n12/dice/rho cannot be trusted or used yet.")
    if log_rows:
        lines.append(
            "- MiXeR run logs contain failure signatures, dominated by: "
            + ", ".join(f"{key}={value}" for key, value in sorted(error_counts.items()))
            + "."
        )
    if missing_current:
        lines.append(
            "- Existing MiXeR sumstats do not cover the current P0 design traits: "
            + ", ".join(missing_current)
            + "."
        )
    lines.append("- The legacy MiXeR scripts target the older gut-only design and do not include ADHD, SCZ, or anxiety benchmarks.")
    lines.append("")
    lines.append("## Sumstats Inventory")
    if sumstats_rows:
        lines.append("| root | trait | size_mb | required_columns | sample_issue |")
        lines.append("|---|---:|---:|---:|---|")
        for row in sumstats_rows:
            issue = row["sample_issue"] or ""
            lines.append(
                f"| {row['root']} | {row['trait']} | {row['size_mb']} | {row['has_required_columns']} | {issue} |"
            )
    else:
        lines.append("- No MiXeR sumstats were found.")
    lines.append("")
    lines.append("## Existing Output Replicate Gate")
    lines.append("- Pass criterion used here: at least 18/20 non-empty replicate JSONs plus a non-empty combined JSON.")
    if expectation_rows:
        failed = [row for row in expectation_rows if row.status == "fail"]
        lines.append(f"- Failed checks: {len(failed)} / {len(expectation_rows)}.")
        for row in failed[:12]:
            lines.append(
                f"- {row.label} {row.phase}: reps={row.nonempty_reps}/{row.expected_reps}, combined={row.combined_present}, root={row.root}"
            )
        if len(failed) > 12:
            lines.append(f"- Additional failed rows omitted from report body: {len(failed) - 12}; see CSV.")
    else:
        lines.append("- No expected-output rows were generated because no configured root exists.")
    lines.append("")
    lines.append("## Log Failures")
    if log_rows:
        for row in log_rows[:20]:
            lines.append(f"- {Path(row['log']).name}:{row['line_number']} {row['line']}")
        if len(log_rows) > 20:
            lines.append(f"- Additional log error rows omitted from report body: {len(log_rows) - 20}; see CSV.")
    else:
        lines.append("- No error signatures were detected in available logs.")
    lines.append("")
    lines.append("## QC Decision")
    lines.append("- Current MiXeR outputs should **not** be used for interpretation or manuscript tables.")
    lines.append("- Required next action: rerun MiXeR from P0-standardized sumstats, then repeat this QC before overlap interpretation.")
    lines.append("- Minimum rerun scope: ASD, IBS, ADHD, SCZ, anxiety, IBD_deLange, CD_deLange, UC_deLange; bivariate pairs should be ASD-centered.")
    lines.append("- Suggested acceptance thresholds: ≥18/20 successful fit and test replicates, finite n1/n2/n12/dice/rho with CI, no Docker/API/EOF errors, and consistency checks against LDSC rg direction.")
    lines.append("")
    lines.append("## Machine-Readable Outputs")
    lines.append(f"- Expected outputs QC: `{OUT_DIR / 'mixer_expected_outputs_qc.csv'}`")
    lines.append(f"- Sumstats inventory: `{OUT_DIR / 'mixer_sumstats_inventory.csv'}`")
    lines.append(f"- Log errors: `{OUT_DIR / 'mixer_run_log_errors.csv'}`")
    lines.append(f"- Combined JSON metrics: `{OUT_DIR / 'mixer_combined_json_metrics.csv'}`")
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    roots = [root for root in [DESKTOP_MIXER_ROOT, T7_MIXER_ROOT] if root.exists()]

    expectation_rows: list[OutputExpectation] = []
    sumstats_rows: list[dict[str, Any]] = []
    log_rows: list[dict[str, Any]] = []
    json_rows: list[dict[str, Any]] = []
    for root in roots:
        expectation_rows.extend(audit_expected_outputs(root))
        sumstats_rows.extend(inventory_sumstats(root))
        log_rows.extend(parse_logs(root))
        json_rows.extend(parse_combined_jsons(root))

    write_csv(OUT_DIR / "mixer_expected_outputs_qc.csv", expectation_rows)
    write_csv(OUT_DIR / "mixer_sumstats_inventory.csv", sumstats_rows)
    write_csv(OUT_DIR / "mixer_run_log_errors.csv", log_rows)
    write_csv(OUT_DIR / "mixer_combined_json_metrics.csv", json_rows)

    report = render_report(expectation_rows, sumstats_rows, log_rows, json_rows)
    report_path = OUT_DIR / "mixer_qc_report.md"
    report_path.write_text(report, encoding="utf-8")
    print(report_path)


if __name__ == "__main__":
    main()
