from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path

from .annotation_mapping import file_sha256, load_label_mapping, mapping_is_approved
from .utils import ensure_dir


def assert_approved(*, report_path: Path, mapping_path: Path, output_path: Path) -> dict:
    if not report_path.is_file():
        raise FileNotFoundError(f"Annotation qualification report not found: {report_path}")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    mapping = load_label_mapping(mapping_path)
    observed_sha256 = file_sha256(mapping_path)
    recorded_sha256 = str(report.get("label_mapping", {}).get("sha256", ""))
    if not bool(report.get("structural_checks_passed", False)):
        raise RuntimeError("Annotation qualification structural checks did not pass.")
    if observed_sha256 != recorded_sha256:
        raise RuntimeError(
            "Annotation mapping changed after the qualification report was generated."
        )
    if not mapping_is_approved(mapping) or not bool(report.get("approved_for_downstream", False)):
        pending = int(mapping["review_status"].ne("approved").sum())
        raise RuntimeError(
            f"Annotation qualification is not approved for downstream analysis: {pending} mapping rows remain."
        )

    payload = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "approved",
        "qualification_report": str(report_path.resolve()),
        "qualification_report_sha256": file_sha256(report_path),
        "mapping_path": str(mapping_path.resolve()),
        "mapping_sha256": observed_sha256,
        "mapping_version": str(mapping["mapping_version"].iloc[0]),
    }
    ensure_dir(output_path.parent)
    temporary = output_path.with_name(output_path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    os.replace(temporary, output_path)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", required=True)
    parser.add_argument("--label-map", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    payload = assert_approved(
        report_path=Path(args.report),
        mapping_path=Path(args.label_map),
        output_path=Path(args.out),
    )
    print(
        f"[assert_annotation_qualification] approved mapping {payload['mapping_version']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
