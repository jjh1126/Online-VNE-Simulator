# Online VNE Simulator

C++17 implementation of VNE-BCP and the VNE-RFD-B, VNE-RFD-D, D-ViNE-LB, and R-ViNE-LB baselines.

## Environment and build

The code was tested on Ubuntu 22.04 with GCC 11.4, CMake 4.3, and GLPK 5.0.

```bash
sudo apt install build-essential cmake libglpk-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Running experiments

Place each case under `examples/<case>/` with matching `substrate_<seed>.txt` and `requests_<seed>.txt` files, then run:

```bash
./build/vne_sim --case <case> --algo VNE-BCP --seed 10 \
  --sample-mode paired --warmup 10000 --measurement-window 50000 \
  --results-dir results
```

`--seed N` loads seeds `0` through `N-1`; use `--seed-offset` to change the first seed. `paired` evaluates matching substrate/request seeds, while `cartesian` evaluates every combination.

Algorithms: `VNE-BCP`, `VNE-RFD-B`, `VNE-RFD-D`, `D-ViNE-LB`, and `R-ViNE-LB`.

Main options: `--threads`, `--horizon`, `--tau`, `--ksp`, `--delta`, `--max-d`, `--alpha`, `--beta`, `--lambda`, `--gamma`, `--vine-timeout`, `--vine-seed`, and `--vine-diagnostics`.

Each run writes `summary.csv`, `runs.csv`, and `requests.csv` to a timestamped result directory.
