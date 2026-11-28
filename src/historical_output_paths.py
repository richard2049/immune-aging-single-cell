"""Keep historical analysis commands from writing into protected outputs."""

from pathlib import Path


def validate_historical_outputs(*paths: str) -> None:
    for path in paths:
        parts = [part.casefold() for part in Path(path).resolve().parts]
        for parent, child in zip(parts, parts[1:]):
            if (parent, child) == ("docs", "assets") or (
                parent == "results" and child.startswith(("blood_age_atlas", "gse164378"))
            ):
                raise ValueError(
                    f"Historical output cannot be written to {path}. "
                    "Choose a separate scratch destination."
                )
