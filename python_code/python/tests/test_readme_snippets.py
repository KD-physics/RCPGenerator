"""Every Python snippet in the README section "Ellipsoids and rods" runs as printed.

    python tests/test_readme_snippets.py            # all snippets of the section (about a minute)
    python tests/test_readme_snippets.py --quick    # the extraction checks and the first snippet
    pytest tests/test_readme_snippets.py            # same; RCP_QUICK=1 for the subset

The section is found by its heading; every ```python block between that heading and the next heading of the same
level is one snippet. Each runs in its own subprocess, in a temporary directory, with a non-interactive Matplotlib
backend and one thread (the pattern of tests/test_examples.py), so each must be complete: imports included, small
particle counts. A snippet fails when it exits with an error or triggers the small-periodic-box warning.

Also run: the complete example in the paragraph of "Quench to force balance" on retrying a quench that stalls.

The other sections of the README hold fragments (`packing = rcpgenerator.Packing(...)`, calls on an existing
`packing`) that are not meant to run alone; the snippets of rcpgenerator.example() are run by tests/test_examples.py.
"""

from __future__ import annotations

import os
import pathlib
import re
import subprocess
import sys
import tempfile
import time

README = pathlib.Path(__file__).resolve().parent.parent / "README.md"
SECTION = "## Ellipsoids and rods"
MIN_SNIPPETS = 10      # the section had 11 when this test was written: fewer means the extraction broke
QUICK = os.environ.get("RCP_QUICK", "") not in ("", "0")

_BLOCK = re.compile(r"^```python[ \t]*\n(.*?)^```[ \t]*$", re.S | re.M)


def section_text() -> str:
    """The README from the section heading up to the next heading of the same level."""
    lines = README.read_text(encoding="utf-8").splitlines(keepends=True)
    starts = [k for k, line in enumerate(lines) if line.rstrip() == SECTION]
    assert len(starts) == 1, f"README must have exactly one heading '{SECTION}', found {len(starts)}"
    end, fenced = len(lines), False
    for k in range(starts[0] + 1, len(lines)):
        if lines[k].startswith("```"):
            fenced = not fenced
        elif not fenced and lines[k].startswith("## "):   # a '#' line inside a code block is a comment
            end = k
            break
    return "".join(lines[starts[0]:end])


def snippets() -> list[tuple[str, str]]:
    """(label, code) for every ```python block of the section; the label names the README line it starts on."""
    text = section_text()
    offset = README.read_text(encoding="utf-8").index(text)
    before = README.read_text(encoding="utf-8")[:offset].count("\n")
    out = []
    for m in _BLOCK.finditer(text):
        line = before + text[:m.start()].count("\n") + 1
        out.append((f"README.md line {line}", m.group(1)))
    return out + retry_snippets()


RETRY_START, RETRY_END = "**If a quench stalls or fails, retry it from the packed state**", "**Where it works.**"


def retry_snippets() -> list[tuple[str, str]]:
    """The complete example of the paragraph on retrying a quench, in the section "Quench to force balance"."""
    full = README.read_text(encoding="utf-8")
    assert full.count(RETRY_START) == 1 and full.count(RETRY_END) == 1, "the quench-retry paragraph of the README was not found"
    a, b = full.index(RETRY_START), full.index(RETRY_END)
    found = [(f"README.md line {full[:a + m.start()].count(chr(10)) + 1} (quench retry)", m.group(1)) for m in _BLOCK.finditer(full[a:b])]
    assert found, "the quench-retry paragraph has no python block"
    return found


def run_snippet(code: str) -> tuple[bool, str, float]:
    env = dict(os.environ, MPLBACKEND="Agg", OMP_NUM_THREADS="1")
    with tempfile.TemporaryDirectory() as tmp:
        t0 = time.time()
        proc = subprocess.run([sys.executable, "-c", code], cwd=tmp, env=env, capture_output=True, text=True,
                              timeout=900)
        return proc.returncode == 0, (proc.stdout + proc.stderr)[-2000:], time.time() - t0


def test_section_has_snippets():
    found = snippets()
    assert len(found) >= MIN_SNIPPETS, f"only {len(found)} python snippets found in '{SECTION}'"
    for label, code in found:
        assert "import rcpgenerator" in code, f"{label}: not self-contained (no 'import rcpgenerator')"
        assert "..." not in code, f"{label}: a snippet with '...' cannot run as printed"
        compile(code, label, "exec")       # a syntax error is reported with the README line
    text = section_text()
    for word in ("aspect_ratios", "orientations", "rod_length", "rod_radius", "rod_capped", "phi_contact",
                 "shrink_to_contact", "ell_permute_orient", "dt_max"):
        assert any(word in code for _, code in found), f"no snippet of the section uses {word}"
        assert word in text


def test_snippets_run():
    found = snippets()
    if QUICK:
        found = found[:1]
    failed = []
    for label, code in found:
        ok, out, dt = run_snippet(code)
        if ok and "periodic box is small" in out:
            ok, out = False, "a snippet should not trigger the small periodic box warning\n" + out
        print(f"  {'ok  ' if ok else 'FAIL'}  {label} ({dt:.1f} s)", flush=True)
        if not ok:
            failed.append(f"{label}:\n{code}\n{out}")
    assert not failed, "\n\n".join(failed)


TESTS = [test_section_has_snippets, test_snippets_run]


def main(argv: list[str]) -> int:
    global QUICK
    QUICK = QUICK or "--quick" in argv
    failures = 0
    for test in TESTS:
        t0 = time.time()
        try:
            test()
            print(f"PASS  {test.__name__}  ({time.time() - t0:.1f} s)", flush=True)
        except Exception as exc:
            failures += 1
            print(f"FAIL  {test.__name__}: {type(exc).__name__}: {exc}", flush=True)
    print(f"{len(TESTS) - failures}/{len(TESTS)} passed" + (" (quick)" if QUICK else ""))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
