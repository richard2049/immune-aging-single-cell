# Troubleshooting (immune_aging_scvi)

## Conda Or Mamba Cannot Reach Channels

Errors such as `HTTP 000`, `Download error (7)`, `WinError 10013`, or failure
to connect to `conda.anaconda.org:443` mean that the environment update did not
run. A subsequent `ModuleNotFoundError` is an expected consequence of the
stale environment, not a separate package-resolution result.

Run these checks in a normal PowerShell terminal:

```powershell
Get-ChildItem Env: | Where-Object { $_.Name -match '^(HTTP_PROXY|HTTPS_PROXY|ALL_PROXY)$' }
Test-NetConnection conda.anaconda.org -Port 443
python -c "import urllib.request as u; r=u.urlopen(u.Request('https://conda.anaconda.org/conda-forge/noarch/repodata.json.zst', method='HEAD'), timeout=15); print(r.status)"
```

If a proxy variable points to a disabled local endpoint such as
`http://127.0.0.1:9`, clear it for the current terminal and retry a read-only
channel query:

```powershell
Remove-Item Env:HTTP_PROXY,Env:HTTPS_PROXY,Env:ALL_PROXY -ErrorAction SilentlyContinue
$mamba = Join-Path ((conda info --base).Trim()) "Library\bin\mamba.exe"
if (-not (Test-Path $mamba)) { throw "mamba.exe not found at $mamba" }
& $mamba repoquery search zlib --override-channels -c conda-forge
```

Do not remove a real institutional proxy. If the TCP test succeeds but the
Python request fails with `WinError 10013`, an application-level firewall,
antivirus, or endpoint-security policy is blocking outbound HTTPS. Allow the
Miniforge `mamba.exe`, `conda.exe`, and `python.exe` executables, or ask the
system administrator to do so. If both HTTPS tests fail, check VPN or network
policy and try another trusted network. `mamba clean -a` is not the first
remedy when every channel is unreachable.

After the read-only query succeeds, run the update with failure guards:

```powershell
$mamba = Join-Path ((conda info --base).Trim()) "Library\bin\mamba.exe"
if (-not (Test-Path $mamba)) { throw "mamba.exe not found at $mamba" }
& $mamba env update -n immune-aging-scvi -f environment.yml --prune
if ($LASTEXITCODE -ne 0) { throw "Environment update failed." }
conda activate immune-aging-scvi
python -c "import scanpy, scvi, celltypist, xgboost; print(scanpy.__version__, scvi.__version__, celltypist.__version__, xgboost.__version__)"
if ($LASTEXITCODE -ne 0) { throw "Required scientific packages are unavailable." }
python -m unittest discover -s tests -v
```

If Bitdefender permits the request only with a broad application rule, use that
rule temporarily for the exact Miniforge executables while updating the
environment. Do not enable "all applications" and do not disable the firewall.
Avoid `Custom Remote Address` unless its address and port semantics have been
verified in the installed Bitdefender version. After package and CellTypist
model downloads complete, remove the temporary broad rules or tighten them one
setting at a time, testing after each change. If selecting TCP causes the rule
not to match, retain `Any` protocol only for the temporary installation rule;
HTTPS still uses TCP, but Bitdefender may classify auxiliary traffic or the
rule differently. Keep the rule application-specific and outbound, never
inbound, and remove it when installation and model downloads finish. Do not
leave a broad outbound rule for a general Python interpreter permanently.

## GPU Training

The tracked configs and `environment.yml` default to CPU for portability and
CI. `environment-gpu.yml` installs the official PyTorch 2.5.1 CUDA 12.1 wheels
in a separate environment; the wheels include their user-space CUDA runtime,
so a separate CUDA Toolkit installation is not required for this workflow.
The NVIDIA display driver must support CUDA 12.1 (Windows driver 527.41 or
newer according to NVIDIA's CUDA 12.1 release notes).

Create and inspect the optional environment:

```powershell
conda env create -n immune-aging-scvi-gpu -f environment-gpu.yml
conda activate immune-aging-scvi-gpu
python -m pip check
python -c "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available()); print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU only')"
nvidia-smi
```

Use `scvi.accelerator: gpu` only in the local config for a qualified run. Use
`scvi.accelerator: cpu` to execute the same CUDA-enabled environment on CPU;
no package reinstall is needed. Do not edit a tracked study profile merely to
switch hardware for a local run.

Before a long run, require a bounded scVI training smoke test on both backends.
Record the PyTorch version, CUDA build, driver, device, elapsed time, and peak
GPU memory. Stop and use CPU if CUDA initialization fails, required DLLs are
missing, memory is exhausted, or the driver resets. Avoid combining Conda CUDA
runtime packages, a CUDA-enabled XGBoost build, and PyTorch wheels built for a
different CUDA release in one environment.

Normal GPU compute is managed by the driver and its thermal/power limits, but
no software environment can guarantee zero hardware wear. Keep vents clear,
use mains power and the manufacturer's normal performance profile, and inspect
temperature and memory with `nvidia-smi` during the first full run. A short
functional smoke test is preferable to an unnecessary prolonged stress test.

### NumPy Or SciPy Linear Algebra Stalls On Windows

Test the numerical backend independently before attributing a stalled
qualification or model to the workflow:

```powershell
python -c "import numpy as np; a=np.array([[1.,2.],[3.,4.]]); print(np.linalg.matrix_rank(a)); print(np.linalg.svd(a, compute_uv=False))"
```

The qualified GPU profile constrains BLAS/LAPACK 3.9 with MKL 2024.2 and its
matching Intel OpenMP and TBB runtime. Recreate or update the GPU environment
from `environment-gpu.yml` if this smoke test hangs. Do not downgrade only MKL
or use `--no-deps`: a partial transaction can leave incompatible OpenMP or TBB
libraries. Preserve an explicit environment export before any repair and
repeat NumPy/SciPy, CPU PyTorch, CUDA, imports, and tests afterwards.

### OpenMP Runtime Conflict On Windows

`OMP Error #15: Initializing libiomp5md.dll, but found libiomp5md.dll already
initialized` means that incompatible OpenMP runtimes were loaded into one
Python process. Do not set `KMP_DUPLICATE_LIB_OK=TRUE`; it suppresses the guard
without making the numerical stack coherent.

Inspect the resolved environment before changing it:

```powershell
conda list -n immune-aging-scvi | Select-String -Pattern "^(mkl|intel-openmp|llvm-openmp|libclang13|numexpr|tbb)\s"
```

The qualified Windows environment uses MKL 2024.2, Intel OpenMP 2024.2, TBB
2021.13, `libclang13 22.1.8 default_hf735972_3`, and `numexpr 2.14.1
mkl_py310h4c4f63a_1`, with no `llvm-openmp` package. Package build strings are
platform-specific, so they are recorded here rather than imposed on the
cross-platform environment file. Preview any repair transaction and reject it
if it removes required scientific packages or upgrades the qualified MKL
stack. After repair, require NumPy/SciPy linear algebra, NumExpr, statsmodels,
CPU PyTorch, bounded CUDA, imports, and the test suite to pass.

### Snakemake Stalls Before Parsing On Windows

If `snakemake --list` or a dry-run stalls before reading the Snakefile, capture
a Python faulthandler trace before rerunning. A trace ending in Snakemake
`SourceCache` while creating a temporary directory under the user-local cache
indicates an endpoint-security or cache-permission problem, not a workflow DAG
failure.

Keep endpoint protection enabled. Permit only the affected user-local
Snakemake cache path, or use a security-policy-approved cache location if the
installed Snakemake version supports it. While that platform issue is being
resolved, run the documented direct modules sequentially from validated
checkpoints, with explicit inputs, outputs, logs, and timeouts. Do not run the
direct modules concurrently with Snakemake and do not use a broad antivirus
exclusion for Python or the repository.

## Windows Cache Permissions

Some Windows endpoint-security configurations allow Python to read the
scientific environment but block temporary files created by Numba while
Scanpy is imported. The visible command can then appear to stall before model
initialization. Keep the security controls enabled and redirect only Numba's
disposable cache to the ignored workspace directory:

```powershell
$numbaCache = Join-Path (Resolve-Path ".tmp").Path "numba-cache"
New-Item -ItemType Directory -Force -Path $numbaCache | Out-Null
$env:NUMBA_CACHE_DIR = $numbaCache
python -c "import scanpy, scvi; print(scanpy.__version__, scvi.__version__)"
```

Set the variable in each new terminal used for Scanpy or scVI. This does not
change model inputs, parameters, or numerical libraries. Do not solve a cache
permission failure by disabling antivirus protection globally or by granting
Python unrestricted filesystem access.

For a stable local environment, the same cache path can be stored as an
environment-scoped Conda variable. Reactivate the environment after setting it:

```powershell
conda env config vars set -n immune-aging-scvi "NUMBA_CACHE_DIR=$numbaCache"
conda deactivate
conda activate immune-aging-scvi
python -c "import celltypist; print(celltypist.__version__)"
```

## Out-of-memory
- Do NOT densify neighbor graphs (`.toarray()`).
- Do NOT store full scVI normalized expression as a dense layer.
- Keep `adata.X` sparse and use layers["counts"] as sparse.

These were common crash points in the original monolithic script.
