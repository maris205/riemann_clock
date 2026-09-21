#!/usr/bin/env python3
"""Run the frozen J0643 main-complex campaign using the shared model engine."""
from pathlib import Path
from archive_complete_model import run_campaign
if __name__=='__main__':
 run_campaign(Path(__file__).resolve().parents[1]/'results/archive_complete/J064326-504112/config.json')
