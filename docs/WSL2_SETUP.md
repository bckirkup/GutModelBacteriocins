# WSL2 CUDA and MPI setup

This setup targets Ubuntu under WSL2 with an NVIDIA GPU. GutIBM also runs as an
MPI/HDF5 CPU build when CUDA is unavailable.

## 1. Update WSL and the Windows NVIDIA driver

From an elevated Windows PowerShell:

```powershell
wsl --update
wsl --shutdown
```

Install or update the Windows NVIDIA driver with WSL2 CUDA support. The Windows
driver is the only display driver required.

**Do not install a Linux NVIDIA display driver inside WSL.** Windows exposes the
CUDA driver to WSL as `libcuda.so`. Installing `cuda`, `cuda-drivers`, or an
ordinary Linux driver meta-package can replace or conflict with that bridge.

NVIDIA guide:
<https://docs.nvidia.com/cuda/wsl-user-guide/>

## 2. Install the WSL-specific CUDA toolkit

Use NVIDIA's CUDA download selector:

<https://developer.nvidia.com/cuda-downloads>

Select:

- Operating System: Linux
- Architecture: x86_64
- Distribution: WSL-Ubuntu
- Version: 2.0
- Installer: deb (network)

Follow the generated commands and install the versioned
`cuda-toolkit-<version>` package. Do not install `cuda`, `cuda-drivers`, or a
Linux display-driver package.

Ensure the toolkit is on `PATH`:

```bash
echo 'export PATH=/usr/local/cuda/bin:$PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=/usr/local/cuda/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}' >> ~/.bashrc
source ~/.bashrc
nvcc --version
nvidia-smi
```

Both commands must work before requesting a CUDA build.

## 3. Install GutIBM build, MPI, HDF5, and Python dependencies

```bash
sudo apt-get update
sudo apt-get install -y \
  build-essential gcc-13 g++-13 cmake \
  libopenmpi-dev openmpi-bin \
  libhdf5-mpi-dev \
  python3 python3-pip python3-venv
```

`rebuild_and_run.sh` creates `.venv` when needed, installs
`python/.[dev]` as an editable package, and uses that interpreter for JSON
discovery and batch runs. To manage it manually:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e "python/.[dev]"
```

If an existing `.venv` has no `pip` module, the helper bootstraps it with
`python -m ensurepip --upgrade` before installing the local package. If that
bootstrap is unavailable, install `python3-venv` and rerun the helper.

If Ubuntu does not provide GCC 13 directly, add the Ubuntu toolchain PPA:

```bash
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install -y gcc-13 g++-13
```

## 4. Build, test, and run

From the repository root:

```bash
./rebuild_and_run.sh
```

The helper:

- caps build parallelism at eight jobs by default to avoid transient WSL2 RAM
  spikes;
- creates or refreshes the repository's `.venv`;
- enables CUDA only when both `nvcc` and a visible GPU are present;
- builds MPI and parallel HDF5 support;
- runs CTest sequentially before experiments;
- launches local MPI ranks with `--bind-to none`.

Useful overrides:

```bash
GUTIBM_BUILD_JOBS=4 ./rebuild_and_run.sh
./rebuild_and_run.sh --cuda off
./rebuild_and_run.sh --cuda on
./rebuild_and_run.sh --mode single \
  --config experiments/smoke_single.json --mpi-ranks 2
./rebuild_and_run.sh --mode batch \
  --config experiments/smoke_batch.json --batch-action dry-run --reuse-build
./rebuild_and_run.sh --mode stage \
  --config experiments/diversity_campaign/stage1_motility_validation \
  --mpi-ranks 1 --reuse-build
./rebuild_and_run.sh --no-gzip-hdf5 --reuse-build --mode single \
  --config experiments/smoke_single.json
```

Successful runs gzip HDF5 outputs to `*.h5.gz` by default (`--no-gzip-hdf5` to
leave them uncompressed).

For a CUDA simulation, its JSON must also set `"gpu_enabled": true`. With more
than one MPI rank and more than one GPU, use `"gpu_device_id": -1` so ranks map
to devices automatically.

## 5. WSL2 memory and CPU limits

WSL2 defaults to 50% of Windows RAM and all logical processors. The reduced GPU
smoke grids should not require increasing that default. If large production
experiments still pressure Windows, create `%UserProfile%\.wslconfig` on the
Windows side and leave enough memory for Windows itself. Example for a 32 GB
host:

```ini
[wsl2]
memory=20GB
processors=8
swap=8GB
```

Apply changes from PowerShell:

```powershell
wsl --shutdown
```

Microsoft reference:
<https://learn.microsoft.com/windows/wsl/wsl-config>

Lower `GUTIBM_BUILD_JOBS` and MPI ranks before increasing the WSL memory cap.
Do not use more MPI ranks than the processors assigned to WSL unless the run is
explicitly an oversubscription test. Stage 1–2 campaign configs are safe at
`mpirun -np 4` after the periodic-x ghost-exchange fix; prefer `--bind-to none`.
**Stage 3** defaults to 1 rank: each rank mirrors the full 50M-cell grid, so
GPU + `np=4` often OOMs WSL and can kill the terminal with no traceback.

## 6. CUDA-aware MPI

Ubuntu's Open MPI packages may not be built with CUDA-buffer support. GutIBM
automatically falls back to host-staged MPI reductions, which is correct but
slower.

Only opt into direct device-buffer collectives when the linked MPI library
reports CUDA support:

```bash
export GUTIBM_CUDA_AWARE_MPI=1
mpirun --bind-to none -np 2 build/tests/test_cuda_aware_mpi_reaction
```

Do not set `GUTIBM_CUDA_AWARE_MPI_FORCE=1` on WSL unless the MPI build is known
to support CUDA buffers.
