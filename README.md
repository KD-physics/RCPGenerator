```markdown
# ND‑RCP: N‑Dimensional Random Close Packing Generator

[![Build Status](https://github.com/<USERNAME>/<REPO>/actions/workflows/build.yml/badge.svg)] [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)]

## Description

ND‑RCP provides fast, flexible tools for generating and relaxing random close packings of spheres in arbitrary dimensions.
- **C++** executables for seeding and optimizing packings with ADAM/Verlet.
- **MATLAB** scripts for initialization, packing, and visualization.

Ideal for computational physicists, ML researchers, and anyone needing controlled particle packings.

## Features

- Seed initial positions at a target packing fraction φ in n dimensions
- Relax packings to high density via gradient‑based dynamics
- Support monodisperse, bidisperse, and power‑law diameter distributions
- Periodic or hard‑walled boundary conditions
- Export/import plain‑text packing files for interoperability
- MATLAB routines for quick prototyping and 2D/3D plotting

## Repository Structure

```

/c++
├ RCPGenerator.cpp         # packing relaxation executable
├ InitializeParticles.cpp  # initial‐positions generator
/matlab
├ initialize\_particlesND.m # N‑D seeding function
├ CreatePacking.m          # packing optimizer
├ plot\_particles\_periodic.m
├ plot\_particles\_3D.m
└ example.m                # end‑to‑end demo
README.md
LICENSE

````

## Prerequisites

- **C++**: GCC or Clang supporting C++17
- **MATLAB**: R2022a or later

## Build & Installation

```bash
cd c++
g++ -O3 -std=c++17 -o InitializeParticles.exe InitializeParticles.cpp
g++ -O3 -std=c++17 -o RCPGenerator.exe       RCPGenerator.cpp
````

*No installation step required for the MATLAB scripts.*

## Usage

### 1. Seed initial positions (C++)

```bash
./InitializeParticles.exe \
  --N 500 \
  --Ndim 3 \
  --phi 0.64 \
  --dist monodisperse \
  --D 1.0 \
  --box 1,1,1 \
  > init_500_3D.txt
```

### 2. Relax packing (C++)

```bash
./RCPGenerator.exe \
  --file init_500_3D.txt \
  --output rcp_500_3D.txt \
  --box 1,1,1 \
  --save-interval 100 \
  --verbose
```

### 3. MATLAB end‑to‑end demo

```matlab
cd matlab
example  % runs initialization, packing, and plots results
```

## Command‑line Arguments (C++)

| Flag              | Type       | Default          | Description                                                     |
| ----------------- | ---------- | ---------------- | --------------------------------------------------------------- |
| `--N`             | int        | —                | Number of particles                                             |
| `--Ndim`          | int        | —                | Number of spatial dimensions                                    |
| `--phi`           | float      | 0.2              | Initial packing fraction                                        |
| `--dist`          | string     | monodisperse     | Diameter distribution: `monodisperse`, `bidisperse`, `powerlaw` |
| `--D`             | float      | 1.0              | Particle diameter (for monodisperse)                            |
| `--D2`            | float      | —                | Second diameter (for bidisperse)                                |
| `--ratio`         | float      | 1.4              | D2/D1 ratio (for bidisperse)                                    |
| `--alpha`         | float      | 3.0              | Exponent (for power‑law distribution)                           |
| `--box`           | comma list | 1,1,…,1          | Box lengths in each dimension (e.g. `--box 1,1,1`)              |
| `--walls`         | bool       | false            | Use hard walls instead of periodic boundaries                   |
| `--save-interval` | int        | 0                | Write intermediate files every N steps (0 = off)                |
| `--verbose`       | bool       | false            | Print progress messages                                         |
| `--file`          | string     | —                | Input file (output of `InitializeParticles.exe`)                |
| `--output`        | string     | `rcp_output.txt` | Output path for the relaxed packing                             |

## MATLAB Functions & Options

### initialize\_particlesND

```matlab
[x, D] = initialize_particlesND(...
    N, phi, box, dist, Name,Value...)
```

| Name      | Type    | Default        | Description                                  |
| --------- | ------- | -------------- | -------------------------------------------- |
| `N`       | int     | —              | Number of particles                          |
| `phi`     | double  | —              | Packing fraction                             |
| `box`     | \[1×n]  | ones(1,n)      | Box lengths in each dimension                |
| `dist`    | string  | 'monodisperse' | 'monodisperse' \| 'bidisperse' \| 'powerlaw' |
| `'D'`     | double  | 1.0            | Diameter (monodisperse)                      |
| `'D2'`    | double  | —              | Second diameter (bidisperse)                 |
| `'ratio'` | double  | 1.4            | D2/D1 ratio (bidisperse)                     |
| `'alpha'` | double  | 3.0            | Exponent (powerlaw)                          |
| `'walls'` | logical | false          | Hard walls if true (otherwise periodic)      |

### CreatePacking

```matlab
[xOpt, DOpt] = CreatePacking(...
    x, D, box, Name,Value...)
```

| Name             | Type   | Default | Description                                       |
| ---------------- | ------ | ------- | ------------------------------------------------- |
| `'maxIter'`      | int    | 10000   | Maximum number of optimization iterations         |
| `'tol'`          | double | 1e-6    | Convergence tolerance                             |
| `'saveInterval'` | int    | 0       | Save intermediate packing every N steps (0 = off) |

## Input / Output Formats

* **Text files**
  Each line:

  ```
  x1 x2 … xn D
  ```

  (space‑separated coordinates + diameter)

* **MATLAB**

  * `x`: N×n matrix of positions
  * `D`: N×1 vector of diameters

## Visualization

* **2D periodic scatter**:

  ```matlab
  plot_particles_periodic(x, D, box)
  ```
* **3D sphere rendering**:

  ```matlab
  plot_particles_3D(x, D, box)
  ```

## Examples

See [`matlab/example.m`](matlab/example.m) for a complete demo:

1. Generate positions
2. Relax packing
3. Plot results

## Contributing

1. Fork the repo
2. Create a feature branch
3. Submit a pull request
4. Address review feedback

## License & Citation

This project is released under the **MIT License**.
If you use ND‑RCP in published work, please cite:

> Desmond, K. *N‑Dimensional Random Close Packing Generator*, GitHub Repository, YYYY.
> DOI or arXiv\:XXXX.XXXXX (if available)

## Contact

Kenneth Desmond — [@<GitHubUsername>](https://github.com/<GitHubUsername>)
[youremail@example.com](mailto:youremail@example.com)

```
```
