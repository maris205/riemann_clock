"""Run the pre-frozen four-line J0530 campaign through the shared engine."""
from pathlib import Path
from archive_complete_model import run_campaign
if __name__ == '__main__':
    base = Path(__file__).resolve().parents[1]
    run_campaign(base / 'results/archive_complete/J053007-250329/selection_frozen.json')
