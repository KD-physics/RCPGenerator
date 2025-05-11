# ND‑RCP: N‑Dimensional Random Close Packing Generator

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

## Description

ND‑RCP provides fast, flexible tools for generating and relaxing random close packings of spheres in arbitrary dimensions.
- **C++** executables for seeding and optimizing packings with ADAM/Verlet.
- **MATLAB** scripts for initialization, packing, and visualization.

Ideal for computational physicists, ML researchers, and anyone needing controlled particle packings.

## Features

- Grows packing diameters in random close packings via ADAM optimizer
- Built-in support for a wide variety of particle size distributions
- Periodic or hard‑walled boundary conditions, including circular or hypersphere walls
- Packings hieight can be constrained to multiple of largest particle diameter
- Export/import plain‑text packing files
- Both c++ and MATLAB routines included

## Repository Structure

/c++

    ├ RCPGenerator.cpp         # packing relaxation executable
    └ InitializeParticles.cpp  # initial‐positions generator

/matlab

    ├ initialize\_particlesND.m # N‑D seeding function
    ├ CreatePacking.m          # packing optimizer
    ├ plot\_particles\_periodic.m
    ├ plot\_particles\_3D.m
    └ example.m                # end‑to‑end demo

README.md

LICENSE


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

## Usage: C++

### 1. Intializing Particle Positions and Diameters

First, particle positions and diameter need to be generated in a column form x, y, z, ..., D. One can do this themselves or used the provide InitializeParticles.cpp function. There is no upper limit on dimensions, only lower limit on 2 dimensions. Note that the positions can be saved to a file and later uploaded using the --file flag in RCPGenerator, or then can be piped into RCPGenerator via | or <. An example use case for InitializeParticles.cpp would be (assuming compiled to InitializeParticles.exe)

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

The command line arguments are as follows

### Command‑line Arguments (C++)

| Flag            | Type       | Default         | Description                                                                                                        |
| --------------- | ---------- | --------------- | ------------------------------------------------------------------------------------------------------------------ |
| `--phi`         | float      | 0.05            | Target packing fraction                                                                                            |
| `--N`           | int        | —               | Number of particles                                                                                                |
| `--Ndim`        | int        | —               | Spatial dimensions (>=2)                                                                                           |
| `--box`         | comma list | 1 repeated Ndim | Box lengths per dimension                                                                                          |
| `--dist`        | string     | —               | Distribution type: mono, bidisperse, gaussian, biGaussian, lognormal, flat, powerlaw, exponential, weibull, custom |
| `--d`           | float      | —               | Diameter (monodisperse)                                                                                            |
| `--d1`          | float      | —               | First diameter (bidisperse)                                                                                        |
| `--d2`          | float      | —               | Second diameter (bidisperse)                                                                                       |
| `--p`           | float      | —               | Fraction (bidisperse, biGaussian)                                                                                  |
| `--mu`          | float      | —               | Mean (gaussian)                                                                                                    |
| `--sigma`       | float      | —               | Std dev (gaussian)                                                                                                 |
| `--mu1`         | float      | —               | Mean1 (biGaussian)                                                                                                 |
| `--sigma1`      | float      | —               | Std1 (biGaussian)                                                                                                  |
| `--mu2`         | float      | —               | Mean2 (biGaussian)                                                                                                 |
| `--sigma2`      | float      | —               | Std2 (biGaussian)                                                                                                  |
| `--d_min`       | float      | —               | Min diameter (flat, powerlaw, exponential)                                                                         |
| `--d_max`       | float      | —               | Max diameter (flat, powerlaw, exponential)                                                                         |
| `--exponent`    | float      | —               | Exponent (powerlaw)                                                                                                |
| `--scale`       | float      | —               | Scale (weibull)                                                                                                    |
| `--shape`       | float      | —               | Shape (weibull)                                                                                                    |
| `--custom_list` | string     | —               | Comma-separated list for custom distribution                                                                       |
| `--fix-height`  | flag       | false           | Fix height dimension when scaling diameters                                                                        |
| `--help`        | flag       | false           | Show help message                                                                                                  |


### 2. Generating RCP (C++)

Using a set of initial positions and diameters provide via column format x,y,z,...D, either piped in or set via the --file flag, and RCP can be generated simply as 

```bash
./RCPGenerator.exe --file init_500_3D.txt
```

where the boundaries are all periodic, and the number of dimensions of particles are inferred from data within the file. With this format the final positions will be printed to the command terminal. If one desires the positions to be printed to a file, the the use case would be

```bash
./RCPGenerator.exe --file init_500_3D.txt --output saved_positions
```

or

```bash
./RCPGenerator.exe --file init_500_3D.txt > saved_positions.txt
```

where .txt is always attached to end of file name automatically. There are many other flags to set the size of the container (box), to make the boundaries hard (walls), to print status updates (verbose), if the final height of the container is to be a multiple of the first listed diameter (fix-height). A full list of options are given below. A more controlled use might be something like 

```bash
./InitializeParticles.exe --N 15000 --Ndim 3 --dist powerlaw --d_min 1.0 --d_max 15.0 --exponent -3 --phi 0.01 --box 1,0.5,1 --walls 0,1,0 > input.txt
./RCPGenerator.exe --file input.txt --output final_packing.txt --NeighborMax 1500 --box 1,0.5,1 --walls 0,1,0 
```

In this example, we have the following attributes
   N (Number of particles)          : 15,000
   Ndim (dimensions)                : 3
   Diameter distribution (--dist)   : Power law with exponent -3 and the lower and upper limits from 1-15
   Packing fraction (phi)           : 0.01
   Container is box of widths (box) : 1 along x, 0.5 along y, and 1 along z
   Boundary conditions (walls)      : 0 (false) means periodic in x, 1 (true) means hard wall in y, 0 (false) means periodic in z
   NeighborMax                      : This one is a little tricky. The code works by building a maintaining a full list of possible nearby neighbors that gets updated periodically. The matrix that stores these possible neighbors is pre-assigned an allocation of memory based on the max number of expected neighbor neighbors. There are default values for NeighborMax but they might not be enough. If it isn't the code will exit and tell you to increase this number. As the size ratio between large and small particles grows, the number of possible nearby neighbors will as well
   
Full list of options

| Flag              | Type       | Default          | Description                                                     |
| ----------------- | ---------- | ---------------- | --------------------------------------------------------------- |
| `--file`          | string     | —                | Input file (output of InitializeParticles.exe)                  |
| `--output`        | string     | packing\_out.txt | Output file for relaxed packing                                 |
| `--box`           | comma list | —                | Box lengths per dimension                                       |
| `--NeighborMax`   | int        | 0 (auto)         | Max neighbors for spatial binning (0 = automatic based on Ndim) |
| `--seed`          | int        | 0                | Seed for RNG (0 = time-based)                                   |
| `--verbose`       | flag       | false            | Print progress and debug messages                               |
| `--fix-height`    | flag       | false            | Fix height dimension when scaling diameters                     |
| `--save-interval` | int        | 0                | Interval (steps) to save intermediate packings (0 = off)        |
| `--walls`         | comma list | 0 repeated Ndim  | Hard-wall flags per dimension (0 = periodic, 1 = hard wall)     |


### 3. How to properly use fixed height


### 4. How to use circular boundary


##. MATLAB

Matlab functionality is one-to-one with c++ code. There is the example file example.m that include end to end demos

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
