# RCPGenerator — Precompiled Python Wheels

**Install RCPGenerator without compiling C++ code.** This folder provides downloadable ZIP archives of precompiled Python wheels (`.whl` files) for the [`rcpgenerator`](https://github.com/KD-physics/RCPGenerator) package.

A wheel is a ready-to-install Python package. **Choose the ZIP for your operating system, extract it, then install the wheel matching your Python version.** You do not need to clone the repository, install a C++ compiler, or build the package from source.

## 1. Check your Python version and computer

Open a terminal (Command Prompt or PowerShell on Windows, Terminal on macOS or Linux) and run:

```bash
python --version
python -c "import platform; print(platform.system(), platform.machine())"
```

If `python` is not recognized on macOS or Linux, try `python3` instead. You need **CPython 3.10, 3.11, 3.12, 3.13, or 3.14** (64-bit, with a supported platform below).

| Your computer | Choose the ZIP/wheels for | Wheel filename contains |
| --- | --- | --- |
| Windows, 64-bit Intel/AMD | Windows x86-64 | `win_amd64` |
| Linux, 64-bit Intel/AMD | Linux x86-64 | `manylinux` and `x86_64` |
| Mac with Apple Silicon (M1/M2/M3/M4, etc.), macOS 15+ | macOS ARM64 / Apple Silicon | `macosx` and `arm64` |
| Mac with Intel processor, macOS 14+ | macOS Intel x86-64 | `macosx` and `x86_64` |

These are the currently documented supported platforms. In particular, a Windows wheel cannot be installed on macOS or Linux, and an Intel Mac wheel is not the same as an Apple Silicon wheel.

## 2. Download the appropriate ZIP

1. In this [`python_precompiled_wheels` folder](https://github.com/KD-physics/RCPGenerator/tree/main/python_precompiled_wheels), click the ZIP archive for your operating system and processor type.
2. On its GitHub file page, choose **Download raw file** (the download icon), and save the ZIP to your computer.
3. **Extract/unzip** the archive (on Windows, right-click → **Extract All**; on macOS, double-click; on Linux, use your archive manager or `unzip`).

You should now have one or more files ending in **`.whl`**. Install a `.whl` file, **not** the downloaded `.zip` archive.

## 3. Select the wheel for your Python version

A wheel filename encodes the Python version and supported platform. Look for the tag that matches the Python version reported in step 1:

| Python version | Wheel tag |
| --- | --- |
| 3.10 | `cp310` |
| 3.11 | `cp311` |
| 3.12 | `cp312` |
| 3.13 | `cp313` |
| 3.14 | `cp314` |

For example, `rcpgenerator-1.0.0-cp312-cp312-win_amd64.whl` is an **illustrative filename** for Python 3.12 on 64-bit Windows. Your actual wheel filename may be different; use the exact filename present in the extracted archive.

## 4. Install the wheel

Open a terminal **in the folder containing your extracted `.whl` files**. Then run:

```bash
python -m pip install "YOUR_MATCHING_WHEEL_FILENAME.whl"
```

Replace `YOUR_MATCHING_WHEEL_FILENAME.whl` with the exact name of your wheel. For example, if the file shown above were present, you would run:

```bash
python -m pip install "rcpgenerator-1.0.0-cp312-cp312-win_amd64.whl"
```

On macOS or Linux, substitute `python3` for `python` if necessary. On Windows, you can also target a particular installed Python version with `py -3.12 -m pip install "YOUR_MATCHING_WHEEL_FILENAME.whl"` (adjust `3.12` to your version).

**That's it.** `pip` installs RCPGenerator into that Python environment, and normally downloads any missing Python dependencies (NumPy, SciPy, pandas, and Matplotlib) automatically. Internet access may be needed for those dependencies. The compiled C++/OpenMP extension is already in the wheel; no local C++ build is required.

> **Using a virtual environment?** Activate it before running the command. **Using Jupyter?** Run `%pip install "/full/path/to/YOUR_MATCHING_WHEEL_FILENAME.whl"` in a notebook cell, ensuring the wheel matches the Python and operating system of the notebook kernel (which may be on a different machine).

## 5. Verify the installation

In the same terminal, run:

```bash
python -c "from rcpgenerator import Packing; print('RCPGenerator installed successfully!')"
```

Or try a small packing in Python:

```python
from rcpgenerator import Packing

p = Packing(
    phi=0.08,
    N=64,
    Ndim=2,
    box=[1.0, 1.0],
    walls=[0, 0],
    dist={"type": "mono", "d": 1.0},
    seed=42,
)
p.pack()
print("Final packing fraction:", p.phi_final)
```

For more examples and the complete Python API, see the [Python package documentation](https://github.com/KD-physics/RCPGenerator/blob/main/python_code/python/README.md) and the [Getting Started notebook](https://github.com/KD-physics/RCPGenerator/blob/main/getting_started.ipynb).

## Troubleshooting

- **`ERROR: ... is not a supported wheel on this platform`** — The wheel does not match your Python version, operating system, or CPU architecture. Recheck `python --version`, the platform command in step 1, and the `cp3xx` / platform tags in the filename.
- **`python` or `pip` not found** — Install a supported Python version from [python.org](https://www.python.org/downloads/). On Windows, try `py -m pip`; on macOS/Linux, try `python3 -m pip`.
- **`No module named rcpgenerator` after installation** — You may have installed into a different Python environment. Use the same interpreter for both `-m pip install` and your script; in notebooks, use `%pip` in the intended kernel.
- **Dependencies cannot be downloaded** — Connect to the internet (or supply those dependencies from your own package source) and retry. The wheel itself contains the compiled RCPGenerator extension, but Python dependencies may still need installation.
- **An old `pip` does not recognize the wheel** — Update it with `python -m pip install --upgrade pip`, then retry.

If there is no wheel for your machine, you can [build from source](https://github.com/KD-physics/RCPGenerator/blob/main/python_code/python/README.md#building-from-source).

---

**Project:** [RCPGenerator](https://github.com/KD-physics/RCPGenerator) · **Questions and bugs:** [GitHub Issues](https://github.com/KD-physics/RCPGenerator/issues)
