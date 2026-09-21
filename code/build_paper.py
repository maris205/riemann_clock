"""Compile the single-column author-review draft with a local TeX installation."""
from pathlib import Path
import os
import subprocess
import shutil

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper"
BUILD = ROOT / "build"
BUILD.mkdir(exist_ok=True)
commands = [
    ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "-output-directory=../build", "main.tex"],
    ["bibtex", "main"],
    ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "-output-directory=../build", "main.tex"],
    ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "-output-directory=../build", "main.tex"],
]
for i, command in enumerate(commands, start=1):
    env = os.environ.copy()
    env["BIBINPUTS"] = str(PAPER) + os.pathsep + env.get("BIBINPUTS", "")
    result = subprocess.run(command, cwd=BUILD if command[0] == "bibtex" else PAPER,
                            env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (BUILD / f"compile_{i}.log").write_text(result.stdout)
    if result.returncode:
        print(result.stdout[-8000:])
        raise SystemExit(result.returncode)
shutil.copyfile(BUILD / "main.pdf", PAPER / "riemann_clock.pdf")
print(PAPER / "riemann_clock.pdf")
