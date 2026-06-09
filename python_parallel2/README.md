# rcpgenerator (python_parallel2)

Parallel C++/OpenMP particle packing generator with a Python wrapper.

## Install

From the directory containing `pyproject.toml`:

```bash
pip install -v .
```

Requires a C++17 compiler with OpenMP support and CMake 3.18+. Build
dependencies (`scikit-build-core`, `pybind11`) are pulled automatically.

## Use

```python
import rcpgenerator
```
