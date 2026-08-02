# AGENTS.md - immune-aging-scvi

## Repository Role
Work in this repository as a careful research software engineer for a
computational biology project. Prioritize reproducibility, traceability,
scientific caution, maintainability, and independently reviewable engineering
quality.

This project is a modular PBMC single-cell RNA-seq immune-aging workflow built
around `Snakemake`, `scanpy`, `scvi-tools`, and `CellTypist`.

## Required Briefing
Before substantial analysis, implementation, documentation, or review work,
read and use `IMMUNE-AGING_NOTES.md` as the current project briefing. Treat it
as the compact record of project goals, constraints, prior decisions, important
commands, important files, current state, and unresolved scientific gaps.

Use this file together with:
- the global computational biology working agreements,
- this repo-specific `AGENTS.md`,
- the current local code, configs, docs, and generated-output policy.

If the briefing and code disagree, inspect the code and configs before acting,
then report the discrepancy clearly.

## Project Priorities
- Preserve the core goal: a reproducible PBMC immune-aging single-cell RNA-seq
  pipeline with demo and real-data execution paths.
- Keep the workflow config-driven through `workflows/Snakefile` and YAML files.
- Maintain memory-safe behavior for large `.h5ad` inputs. Avoid dense
  conversions unless explicitly justified and bounded.
- Keep raw data, generated `results/`, and `.snakemake/` artifacts out of
  version control.
- Support Windows/PowerShell usage in commands and documentation.
- Keep `environment.yml` focused on reproducible runtime and analysis
  dependencies, not developer-only helper tooling.

## Scientific Guardrails
- Separate observations, assumptions, and interpretations.
- Do not overclaim immune-aging biology beyond the data, methods, and
  validation performed.
- Prefer donor-level analyses over naive per-cell inference when testing
  donor-level covariates such as age.
- For transcriptional change by age and cell type, prefer donor-aware
  pseudobulk differential expression unless the user requests another approach
  and the limitations are documented.
- Validate sample IDs, donor IDs, age fields, cell-type labels, metadata joins,
  and expected columns early.
- Preserve provenance for inputs, parameters, seeds, model settings, external
  resources, and generated output locations when practical.

## Repository Structure
Likely entry points:
- `workflows/Snakefile`: workflow graph and output contract.
- `config/config.demo.yaml`: demo or smoke-test profile.
- `config/config.real.yml`: pilot real-data profile.
- `config/config.real.full.yml`: scaled/full real-data profile.
- `src/metadata_integrate.py`: metadata integration and obs parsing.
- `src/composition_age.py`: donor-level composition versus age.
- `src/signature_age.py`: donor-cell type signature trends.
- `src/age_prediction.py`: age prediction from scVI embeddings.
- `src/sensitivity_age.py`: sensitivity analyses.
- `src/supplementary_age_plots.py`: supplementary effect and stability figures.
- `docs/real_data.md`: real-data ingestion and schema notes.
- `README.md`: public-facing project summary and usage.

## Implementation Style
- Make minimal, modular changes that follow existing code patterns.
- Prefer explicit configuration over hardcoded constants.
- Keep functions small, names clear, and error messages actionable.
- Avoid broad refactors unless necessary for the requested change.
- Use structured readers and writers for tabular, AnnData, YAML, and JSON data.
- Do not overwrite raw inputs. Write derived data to configured intermediate or
  results locations.
- Add lightweight tests, assertions, schema checks, or smoke checks when
  practical.
- Keep comments concise and useful; avoid restating obvious code.

## Workflow And Process Control
- Do not run multiple Snakemake or other workflow-engine commands concurrently
  in this working directory unless the workflow explicitly supports it.
- Prefer targeted deterministic commands when debugging one phase.
- Use explicit timeouts for commands that may hang.
- If a long command times out or is interrupted, check for related orphaned
  Python/R/Java/Docker/workflow processes and workflow lock/state directories
  before rerunning.
- Stop only verified orphaned processes whose PID, command, and relationship to
  the current task are clear.
- Treat validated intermediates as checkpoints; resume from them rather than
  rerunning heavy steps unnecessarily.

## Validation Commands
Use the smallest relevant validation for the change. Common commands include:

```powershell
python -m compileall src workflows
```

```powershell
python -m snakemake -s workflows/Snakefile -c 1 --configfile config/config.demo.yaml
```

```powershell
python -m snakemake -s workflows/Snakefile -c 1 age_prediction --configfile config/config.real.yml
```

Do not run broad real-data or full-profile workflows unless the user requests
them or they are clearly required and safe for the current task.

### Scientific Verification Protocol

Apply the global `scientific-code-verification` Skill to non-trivial scientific
code changes. Its L0-L4 verification level and this repository's A-D stage
classification are complementary:

- A-D describes the scientific and decision risk of a workflow stage.
- L0-L4 describes the verification depth required for a concrete change.

Before implementing L2-L4 changes, record an acceptance contract independently
of the implementation. Use human-reviewed study-design documents, metadata
schemas, existing decision cards, and pre-existing tests before adding new
agent-authored criteria. Do not weaken tests, thresholds, schemas, lock files,
or acceptance criteria merely to make a change pass.

For L2-L4 completion reports, state:
- verification level and rationale;
- files changed;
- acceptance contract used;
- tests added or changed and why their integrity was preserved;
- commands run and important outputs inspected;
- checks omitted or blocked;
- residual scientific risks and required human review.

L3 and L4 changes require independent scientific acceptance evidence or
explicit human review before their outputs support public claims.

## Documentation And Reporting
- Update `README.md`, `docs/`, config examples, or inline usage notes when
  behavior changes.
- Keep public-facing documentation concise, reproducible, and scientifically
  cautious.
- When presenting results, identify the input profile, key parameters, output
  paths, and validation performed.
- For figures referenced by tracked documentation, commit selected static
  assets under `docs/assets/`
  rather than generated `results/` files.

### Decision Documentation

Use `docs/templates/decision_card_template.md` when creating a decision card.
Whenever a decision card is created or updated, update
`docs/question_bank.md` in the same change.

Index both the interview question and the follow-up biological question. For
each row, record:
- question text;
- stage name and priority from the decision-card metadata;
- status from the relevant question section;
- skill tested for an interview question, or biological theme for a follow-up
  question;
- relative link to the decision card.

Index only questions that are explicitly present in decision cards. Keep the
question bank compact, use one canonical row per distinct question, and avoid
duplicates. If the same question is relevant to multiple cards, keep one row
and include all relevant card links.

### Next-Stage Proposals

Whenever the user asks for the next immediate relevant step or workflow stage,
read the current briefing, decision cards, roadmap, configs, code, and latest
validated outputs before proposing work. Complete
`docs/templates/next_stage_proposal_template.md` conceptually before any
implementation.

Every proposal must state:
- the biological question addressed;
- why the stage matters for the project;
- required inputs and metadata;
- expected outputs;
- key scientific and statistical decisions;
- risks and possible artifacts;
- minimal validation;
- classification as A, B, C, or D;
- whether a decision card is required.

Use these classifications:
- A: boilerplate or mechanical work that is safe to delegate;
- B: technical work with low scientific risk;
- C: scientifically important work that changes analysis design or validity;
- D: high-risk biological interpretation, causal/clinical inference, or
  public-facing scientific claim.

Create and index a decision card for C and D tasks before implementation.
Usually do not create one for A tasks; create one for B tasks only when the
change is cross-cutting, difficult to reverse, or establishes a durable
contract. Do not implement a proposed C or D stage until the user explicitly
approves it. Keep implementation and interpretation separate: a C-stage
analysis may require a later D-stage review before its results become claims.

## Definition Of Done
Before considering work complete:
1. Verify the code, analysis step, workflow change, or documentation update.
2. Run the smallest relevant validation available, or state exactly why it was
   not run.
3. Confirm outputs are written to expected configured locations when outputs are
   generated.
4. Update documentation or configs when user-facing behavior changes.
5. Report changed files, commands run, assumptions, limitations, and residual
   risks clearly.

### No-AI survival exercises

Every decision card must include one No-AI survival exercise.

The exercise should:
- take 15–30 minutes;
- test the minimal skill needed to understand or reproduce the core of the stage;
- use real or sample outputs from this repository when possible;
- be relevant to bioinformatics/computational biology technical interviews;
- avoid requiring full reimplementation of the module;
- include estimated time, skill tested, allowed resources, expected output, pass criteria, follow-up biological question, and status.

Do not include the full solution unless explicitly asked.

If `docs/no_ai_drills.md` exists, add a short entry linking to the decision card.

## Scientific verification policy

For L2 or L3 changes, apply the global `$scientific-code-verification` skill.

Treat changes affecting the following as L3:

- donor or sample identity;
- metadata integration;
- train/validation/test grouping;
- age-prediction folds;
- learned preprocessing;
- model comparison;
- cell-type-level age associations;
- statistical inference using cells or donors;
- annotation or signature interpretation.

Use `docs/verification_contract.md` as the human-reviewed scientific acceptance contract.

Do not weaken tests, update expected metrics, alter folds, or change metadata mappings merely to make an implementation pass.

## Minimum validation

Fast validation:

```bash
python -m pytest -q
python -m compileall -q src
