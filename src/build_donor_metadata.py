from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd

from .utils import ensure_dir


def _norm_text(x: object) -> str:
    return re.sub(r"[^A-Z0-9]+", "", str(x).upper())


def _first_non_null(values: pd.Series):
    vals = values.dropna()
    if vals.empty:
        return np.nan
    return vals.iloc[0]


def _parse_age(values: pd.Series) -> pd.Series:
    age = values.astype(str).str.extract(r"(\d+(?:\.\d+)?)", expand=False)
    return pd.to_numeric(age, errors="coerce")


def _parse_sex(values: pd.Series) -> pd.Series:
    x = values.astype(str).str.strip().str.lower()
    out = pd.Series(np.nan, index=values.index, dtype=object)
    out[x.str.startswith("m")] = "male"
    out[x.str.startswith("f")] = "female"
    return out


def _extract_sample_and_donor(
    obs_names: pd.Index,
    sample_delim: str,
    donor_index: int,
) -> pd.DataFrame:
    s = pd.Series(obs_names.astype(str), index=obs_names)
    sample = s.str.split("_", n=1).str[0]
    if sample_delim.lower() == "auto":
        candidates = ["-", ".", "|", ":"]
        best_delim = "-"
        best_hits = -1
        probe = sample.head(min(100000, sample.shape[0]))
        for d in candidates:
            hits = int(probe.str.contains(re.escape(d), regex=True).sum())
            if hits > best_hits:
                best_hits = hits
                best_delim = d
        sample_delim = best_delim

    donor = sample.str.split(sample_delim, regex=False).str[donor_index]
    return pd.DataFrame(
        {
            "cell_id": s.to_numpy(),
            "sample_id": sample.to_numpy(),
            "donor_id": donor.to_numpy(),
            "batch": sample.to_numpy(),
        }
    )


def _aggregate_base_metadata(sample_df: pd.DataFrame) -> pd.DataFrame:
    sample_df = sample_df.dropna(subset=["donor_id"]).copy()
    sample_df["donor_id"] = sample_df["donor_id"].astype(str)

    n_cells = sample_df["donor_id"].value_counts().rename("n_cells")
    n_samples = sample_df.groupby("donor_id")["sample_id"].nunique().rename("n_samples")
    batches = (
        sample_df.groupby("donor_id")["batch"]
        .apply(lambda x: ";".join(sorted({str(v) for v in x if pd.notna(v)})))
        .rename("batches")
    )
    out = pd.concat([n_cells, n_samples, batches], axis=1).reset_index(names="donor_id")
    return out.sort_values("donor_id").reset_index(drop=True)


def _read_table(path: Path, sep: str | None, sheet: str | None) -> pd.DataFrame:
    suffix = "".join(path.suffixes).lower()
    if suffix.endswith((".xlsx", ".xls")):
        return pd.read_excel(path, sheet_name=sheet or 0)
    if sep is not None:
        return pd.read_csv(path, sep=sep)
    if suffix.endswith(".tsv") or suffix.endswith(".tsv.gz") or suffix.endswith(".txt"):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def _infer_named_col(df: pd.DataFrame, aliases: Iterable[str]) -> str | None:
    aliases = {_norm_text(a) for a in aliases}
    cols = list(df.columns)
    norm_map = {c: _norm_text(c) for c in cols}

    for c in cols:
        if norm_map[c] in aliases:
            return c
    for c in cols:
        if any(a in norm_map[c] for a in aliases):
            return c
    return None


def _infer_overlap_col(df: pd.DataFrame, target_tokens: set[str]) -> str | None:
    best_col = None
    best_overlap = 0
    for c in df.columns:
        if not (pd.api.types.is_object_dtype(df[c]) or pd.api.types.is_string_dtype(df[c])):
            continue
        vals = df[c].dropna().astype(str)
        if vals.empty:
            continue
        tokens = {_norm_text(v) for v in vals.unique()}
        overlap = len(tokens & target_tokens)
        if overlap > best_overlap:
            best_overlap = overlap
            best_col = c
    return best_col


def _build_partial_from_table(
    df: pd.DataFrame,
    donor_ids: set[str],
    sample_delim: str,
    donor_index: int,
    donor_col: str | None,
    sample_col: str | None,
) -> pd.DataFrame | None:
    df = df.copy()
    if donor_col is None and sample_col is None:
        donor_col = _infer_named_col(
            df,
            aliases=[
                "donor_id",
                "donor",
                "subject",
                "subject_id",
                "participant",
                "participant_id",
                "individual",
            ],
        )
    if sample_col is None:
        sample_col = _infer_named_col(df, aliases=["sample_id", "sample", "batch", "pool", "tube"])

    norm_donor_ids = {_norm_text(x) for x in donor_ids}

    if donor_col is None:
        donor_col = _infer_overlap_col(df, norm_donor_ids)
    if donor_col is not None:
        donor = df[donor_col].astype(str).str.strip()
    elif sample_col is not None:
        donor = df[sample_col].astype(str).str.split(sample_delim, regex=False).str[donor_index]
    else:
        return None

    out = pd.DataFrame({"donor_id": donor.astype(str)})
    out = out[out["donor_id"].str.len() > 0].copy()

    age_col = _infer_named_col(df, aliases=["age", "donor_age", "age_years", "age_yrs", "chronological_age"])
    sex_col = _infer_named_col(df, aliases=["sex", "gender"])
    cohort_col = _infer_named_col(df, aliases=["cohort", "group", "condition", "status"])

    if age_col is not None:
        out["age"] = _parse_age(df.loc[out.index, age_col])
    if sex_col is not None:
        out["sex"] = _parse_sex(df.loc[out.index, sex_col])
    if cohort_col is not None:
        out["cohort"] = df.loc[out.index, cohort_col].astype(str)
    if sample_col is not None:
        out["sample_id_meta"] = df.loc[out.index, sample_col].astype(str)

    keep = [c for c in ["donor_id", "age", "sex", "cohort"] if c in out.columns]
    if len(keep) == 1:
        return None
    out = out[keep].copy()
    out = out.replace({"": np.nan})
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True, help="Input AnnData used by the pipeline.")
    ap.add_argument("--out", required=True, help="Output donor metadata CSV.")
    ap.add_argument(
        "--supp",
        action="append",
        default=[],
        help="Supplementary metadata file (repeatable): csv/tsv/txt/xlsx.",
    )
    ap.add_argument("--supp-sep", default=None, help="Delimiter override for supplementary csv/tsv.")
    ap.add_argument("--supp-sheet", default=None, help="Excel sheet name/index for supplementary xlsx.")
    ap.add_argument(
        "--sample-delim",
        default="auto",
        help="Delimiter in sample token (default: auto; tries '-', '.', '|', ':').",
    )
    ap.add_argument("--donor-index", type=int, default=1, help="Index of donor token after split by sample-delim.")
    ap.add_argument("--donor-col", default=None, help="Explicit donor column name in supplementary files.")
    ap.add_argument("--sample-col", default=None, help="Explicit sample column name in supplementary files.")
    args = ap.parse_args()

    h5ad_path = Path(args.h5ad)
    out_path = Path(args.out)
    ensure_dir(out_path.parent)

    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        sample_df = _extract_sample_and_donor(
            adata.obs_names,
            sample_delim=str(args.sample_delim),
            donor_index=int(args.donor_index),
        )
    finally:
        if hasattr(adata, "file") and getattr(adata, "file", None) is not None:
            adata.file.close()
    donor_meta = _aggregate_base_metadata(sample_df)

    donor_ids = set(donor_meta["donor_id"].astype(str).unique())
    partials: list[pd.DataFrame] = []

    for supp in args.supp:
        path = Path(supp)
        if not path.exists():
            raise FileNotFoundError(f"Supplementary file not found: {path}")
        table = _read_table(path, sep=args.supp_sep, sheet=args.supp_sheet)
        partial = _build_partial_from_table(
            table,
            donor_ids=donor_ids,
            sample_delim=str(args.sample_delim),
            donor_index=int(args.donor_index),
            donor_col=args.donor_col,
            sample_col=args.sample_col,
        )
        if partial is not None and not partial.empty:
            partials.append(partial)

    if partials:
        extra = pd.concat(partials, axis=0, ignore_index=True)
        extra = (
            extra.groupby("donor_id", as_index=False)
            .agg({c: _first_non_null for c in extra.columns if c != "donor_id"})
            .reset_index(drop=True)
        )
        donor_meta = donor_meta.merge(extra, on="donor_id", how="left")

    donor_meta["has_age"] = donor_meta["age"].notna() if "age" in donor_meta.columns else False
    donor_meta["has_sex"] = donor_meta["sex"].notna() if "sex" in donor_meta.columns else False
    donor_meta.to_csv(out_path, index=False)

    print(f"Wrote donor metadata: {out_path}")
    print(f"Donors: {donor_meta.shape[0]}")
    if "age" in donor_meta.columns:
        print(f"Donors with age: {int(donor_meta['age'].notna().sum())}")


if __name__ == "__main__":
    main()
