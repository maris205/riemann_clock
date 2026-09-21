"""Compile the focused, single-column spectroscopy author-review manuscript."""
from pathlib import Path
import os
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / 'paper_spectroscopy'
BUILD = ROOT / 'build_spectroscopy'
BUILD.mkdir(exist_ok=True)
commands = [
    ['pdflatex', '-interaction=nonstopmode', '-halt-on-error', '-output-directory=../build_spectroscopy', 'main.tex'],
    ['bibtex', 'main'],
    ['pdflatex', '-interaction=nonstopmode', '-halt-on-error', '-output-directory=../build_spectroscopy', 'main.tex'],
    ['pdflatex', '-interaction=nonstopmode', '-halt-on-error', '-output-directory=../build_spectroscopy', 'main.tex'],
    ['pdflatex', '-interaction=nonstopmode', '-halt-on-error', '-output-directory=../build_spectroscopy', 'main.tex'],
]
for i, command in enumerate(commands, 1):
    env = os.environ.copy()
    env['BIBINPUTS'] = str(PAPER) + os.pathsep + env.get('BIBINPUTS', '')
    result = subprocess.run(command, cwd=BUILD if command[0] == 'bibtex' else PAPER,
                            env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (BUILD / f'compile_{i}.log').write_text(result.stdout)
    if result.returncode:
        print(result.stdout[-9000:])
        raise SystemExit(result.returncode)
shutil.copyfile(BUILD / 'main.pdf', PAPER / 'feii_relative_frequencies.pdf')
print(PAPER / 'feii_relative_frequencies.pdf')
