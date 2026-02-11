# Troubleshooting (immune_aging_scvi)

## GPU training fails
If `scvi_train.py` fails with CUDA errors, set:
- `config.yaml -> scvi -> accelerator: cpu`

Install a CUDA-enabled torch build only if you know your CUDA version matches.

## Out-of-memory
- Do NOT densify neighbor graphs (`.toarray()`).
- Do NOT store full scVI normalized expression as a dense layer.
- Keep `adata.X` sparse and use layers["counts"] as sparse.

These were common crash points in the original monolithic script.