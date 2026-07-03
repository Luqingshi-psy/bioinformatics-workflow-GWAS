import csv
import gzip
import json
import math
from pathlib import Path


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
AMBIGUOUS = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def safe_float(value):
    if value is None:
        return math.nan
    value = str(value).strip()
    if value in {"", ".", "NA", "NaN", "nan", "null"}:
        return math.nan
    try:
        return float(value)
    except ValueError:
        return math.nan


def finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def complement(allele):
    return allele.translate(COMPLEMENT).upper()


def is_snp(a1, a2):
    return len(a1) == 1 and len(a2) == 1 and a1 in "ACGT" and a2 in "ACGT"


def is_ambiguous(a1, a2):
    return (a1.upper(), a2.upper()) in AMBIGUOUS


def parse_vcf_sample(format_keys, sample_value):
    keys = format_keys.split(":")
    values = sample_value.split(":")
    return dict(zip(keys, values))


def p_from_lp(lp):
    lp = safe_float(lp)
    if not finite(lp):
        return math.nan
    if lp < 0:
        return math.nan
    if lp > 323:
        return 0.0
    return 10 ** (-lp)


def z_from_beta_se(beta, se):
    beta = safe_float(beta)
    se = safe_float(se)
    if not finite(beta) or not finite(se) or se <= 0:
        return math.nan
    return beta / se


def write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def dict_writer(path, fieldnames, gz=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if gz:
        handle = gzip.open(path, "wt", newline="")
    else:
        handle = path.open("w", newline="")
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    return handle, writer

