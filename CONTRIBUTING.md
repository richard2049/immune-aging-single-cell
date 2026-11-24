# Contributing

Contributions that improve reproducibility, validation, portability, or the
scientific clarity of this workflow are welcome. Changes should remain focused
and should preserve the distinction between exploratory analyses and supported
scientific conclusions.

## Development setup

Create the scientific runtime from the tracked environment:

```powershell
conda env create -f environment.yml
conda activate immune-aging-scvi
```

Install development-only tools in a separate virtual environment:

```powershell
$basePython = Join-Path ((conda info --base).Trim()) "python.exe"
& $basePython -m venv .venv-dev
$devPython = Join-Path (Resolve-Path ".venv-dev") "Scripts\python.exe"
& $devPython -m pip install -r requirements-dev.txt
```

## Development principles

- Keep workflow behavior configuration-driven and avoid hidden local state.
- Treat raw inputs and validated intermediates as immutable.
- Keep generated data, results, models, environments, and caches out of Git.
- Validate sample, donor, replicate, age, and cell-type metadata explicitly.
- Use donor-aware inference for donor-level variables such as age.
- Document assumptions and limitations without converting associations into
  causal, longitudinal, or clinical claims.

## Validation

Run the smallest checks relevant to the change. The standard fast checks are:

```powershell
python -m unittest discover -s tests -v
python -m compileall -q src
& $devPython -m ruff check src tests workflows
& $devPython -m ruff format --check src tests workflows
```

Workflow changes should also construct the relevant Snakemake DAG before a
full execution. Scientific changes must follow the acceptance criteria in the
[scientific verification contract](docs/verification_contract.md). Do not
weaken tests, schemas, grouping rules, or expected outputs merely to make a
change pass.

## Pull requests

A pull request should describe its scope, the validation performed, and any
remaining scientific or technical limitations.

Use the public repository for code, documented methods and provenance,
environment specifications, configuration examples, tests, and reviewed
results that help readers reproduce or assess the workflow. Cite the accession
and document how to obtain the inputs under their access terms. Store private
notes, intermediate reviews, participant-level extracts, generated results,
and local configuration under ignored paths such as `.local/`, `data/`, or
`results/`.

Review the paths and contents of every outgoing commit before pushing to any
public branch or tag. Check PR text, attachments, and release files before
publishing. Publish a figure or table only after reviewing its provenance,
data permissions, identifiers, and scientific claims. Add new documentation
assets to the ignore exceptions by exact path.
