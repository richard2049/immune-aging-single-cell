from __future__ import annotations

import argparse
import os
import re
import warnings
from pathlib import Path
from typing import Iterable

from .utils import ensure_dir

try:
    import synapseclient
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "synapseclient is not installed in this environment. "
        "Install it with `pip install synapseclient` in your active env."
    ) from exc


FOLDER_TYPES = {
    "org.sagebionetworks.repo.model.Project",
    "org.sagebionetworks.repo.model.Folder",
}
FILE_TYPES = {
    "org.sagebionetworks.repo.model.FileEntity",
    "org.sagebionetworks.repo.model.table.MaterializedView",
}


def _split_csv_values(text: str | None) -> list[str]:
    if not text:
        return []
    return [x.strip().lower() for x in text.split(",") if x.strip()]


def _login(profile: str | None, auth_token_env: str) -> synapseclient.Synapse:
    syn = synapseclient.Synapse()
    token = os.getenv(auth_token_env, "").strip()
    if token:
        syn.login(authToken=token, silent=True)
        return syn
    if profile:
        syn.login(profile=profile, silent=True)
        return syn
    syn.login(silent=True)
    return syn


def _iter_children(
    syn: synapseclient.Synapse,
    root_id: str,
    recursive: bool,
) -> Iterable[dict]:
    stack = [root_id]
    while stack:
        current = stack.pop()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            for child in syn.getChildren(current):
                yield child
                if recursive and child.get("type") in FOLDER_TYPES:
                    stack.append(child["id"])


def _safe_name(text: str) -> str:
    text = re.sub(r"[^\w.\-]+", "_", text).strip("_")
    return text or "entity"


def _matches_keywords(
    name: str,
    include_keywords: list[str],
    exclude_keywords: list[str],
) -> bool:
    lname = name.lower()
    if include_keywords and not any(k in lname for k in include_keywords):
        return False
    if exclude_keywords and any(k in lname for k in exclude_keywords):
        return False
    return True


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--entity-id", default="syn56693935", help="Synapse project/folder/entity ID.")
    ap.add_argument("--recursive", action="store_true", help="Traverse subfolders recursively.")
    ap.add_argument(
        "--keywords",
        default="meta,metadata,clinical,donor,age,phenotype,covariate",
        help="Comma-separated include keywords to match in child names.",
    )
    ap.add_argument(
        "--exclude",
        default="raw,counts,matrix,h5ad,mtx,loom,bam,fastq",
        help="Comma-separated exclude keywords to avoid large expression files.",
    )
    ap.add_argument("--profile", default=None, help="Synapse profile name in .synapseConfig.")
    ap.add_argument(
        "--auth-token-env", default="SYNAPSE_AUTH_TOKEN", help="Env var name for Synapse PAT."
    )
    ap.add_argument(
        "--download-dir", default="data/raw/raw_counts_h5ad", help="Destination folder."
    )
    ap.add_argument("--download", action="store_true", help="Download matched files.")
    ap.add_argument(
        "--limit", type=int, default=0, help="Max number of matched files (0 = no limit)."
    )
    args = ap.parse_args()

    include_keywords = _split_csv_values(args.keywords)
    exclude_keywords = _split_csv_values(args.exclude)

    syn = _login(profile=args.profile, auth_token_env=args.auth_token_env)
    ensure_dir(args.download_dir)

    matches: list[dict] = []
    for child in _iter_children(syn, args.entity_id, recursive=bool(args.recursive)):
        ctype = child.get("type", "")
        if ctype not in FILE_TYPES:
            continue
        name = str(child.get("name", ""))
        if not _matches_keywords(name, include_keywords, exclude_keywords):
            continue
        matches.append(child)
        if args.limit and len(matches) >= args.limit:
            break

    if not matches:
        print("No matching file entities found with current filters.")
        print(f"Entity: {args.entity_id}")
        print(f"Include keywords: {include_keywords}")
        print(f"Exclude keywords: {exclude_keywords}")
        return

    print(f"Matched {len(matches)} candidate metadata files:")
    for idx, m in enumerate(matches, start=1):
        print(f"{idx:3d}. {m['id']}  {m.get('name', '')}")

    if not args.download:
        print("\nDry run only. Add --download to fetch matched files.")
        return

    dest = Path(args.download_dir)
    for m in matches:
        syn_id = m["id"]
        name = _safe_name(m.get("name", syn_id))
        print(f"Downloading {syn_id} ({name})...")
        entity = syn.get(syn_id, downloadLocation=str(dest))
        local_path = getattr(entity, "path", None)
        if local_path:
            print(f"  -> {local_path}")


if __name__ == "__main__":
    main()
