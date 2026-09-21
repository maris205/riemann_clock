#!/usr/bin/env python3
"""Archive small NIST ASD queries and parse an auditable transition lookup.

Default: use existing snapshots (no network if all snapshots exist).
--refresh: explicitly retrieve current NIST output; the remote database may change.
Only the Python standard library is required. No downloaded code is executed.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
import urllib.parse
import urllib.request

BASE = Path(__file__).resolve().parents[1]
DEST = BASE / "data" / "atomic"


def query(ion: str, lo: int, hi: int, fmt: str = "3") -> str:
    p = dict(spectra=ion, low_w=str(lo), upp_w=str(hi), unit="0", format=fmt,
             line_out="0", en_unit="0", output="0", page_size="100",
             show_obs_wl="1", show_calc_wl="1", unc_out="1", order_out="0",
             show_av="3", tsb_value="0", A_out="0", intens_out="on",
             enrg_out="on", conf_out="on", term_out="on", J_out="on",
             allowed_out="1", forbid_out="1", max_low_enrg="1", min_accur="",
             bibrefs="1")
    return "https://physics.nist.gov/cgi-bin/ASD/lines1.pl?" + urllib.parse.urlencode(p)


SOURCES = {
    "nist_fe_ii_ground_lines.tsv": query("Fe II", 1500, 2700),
    "nist_mg_ii_ground_lines.tsv": query("Mg II", 2795, 2805),
    "nist_fe_ii_query.html": query("Fe II", 1500, 2700, "0"),
    "nist_asd_lineshelp.html": "https://physics.nist.gov/PhysRefData/ASD/Html/lineshelp.html",
}


def table(filename: str) -> list[dict]:
    rows = list(csv.DictReader(io.StringIO((DEST / filename).read_text()), delimiter="\t"))
    return [{str(k): str(v).strip() for k, v in row.items() if k} for row in rows
            if row.get("ritz_wl_vac(A)")]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args()
    DEST.mkdir(parents=True, exist_ok=True)
    manifest_path = DEST / "source_manifest.json"
    old = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    old_entries = {r["file"]: r for r in old.get("sources", [])}
    entries = []
    for filename, url in SOURCES.items():
        path = DEST / filename
        refreshed = args.refresh or not path.exists()
        if refreshed:
            request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (research data archival)"})
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = response.read(5_000_001)
            if len(payload) > 5_000_000:
                raise RuntimeError(f"Unexpectedly large NIST response: {filename}")
            if filename.endswith("tsv") and not payload.startswith(b"obs_wl_vac(A)"):
                raise RuntimeError(f"Not a vacuum-wavelength NIST table: {filename}")
            path.write_bytes(payload)
        payload = path.read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        previous = old_entries.get(str(path.relative_to(BASE)), {})
        if not refreshed and previous.get("sha256", digest) != digest:
            raise RuntimeError(f"Archived source checksum changed: {filename}")
        entries.append({
            "file": str(path.relative_to(BASE)), "url": url, "bytes": len(payload),
            "sha256": digest,
            "retrieved_utc": datetime.now(timezone.utc).isoformat() if refreshed else
                previous.get("retrieved_utc", "2026-09-20; initial retrieval before manifest creation"),
            "access": "Publicly accessible official NIST ASD; no authentication",
        })
    html = (DEST / "nist_fe_ii_query.html").read_text()
    version_match = re.search(r"ver\.&nbsp;([0-9.]+)", html)
    raw_fe = table("nist_fe_ii_ground_lines.tsv")
    raw_mg = table("nist_mg_ii_ground_lines.tsv")
    selected = []
    for ion, data, targets in [("Fe II", raw_fe, [1608, 1611, 2249, 2260, 2344, 2374, 2382, 2586, 2600]),
                               ("Mg II", raw_mg, [2796, 2803])]:
        for target in targets:
            matches = [r for r in data if r.get("Type", "") == ""
                       and abs(float(r["ritz_wl_vac(A)"]) - target) < 1.5]
            if len(matches) != 1:
                raise RuntimeError(f"Non-unique allowed transition {ion} {target}: {len(matches)}")
            r = matches[0]
            wl = float(r["ritz_wl_vac(A)"])
            err = r.get("unc_ritz_wl", "")
            assert float(r["Ei(cm-1)"]) == 0.0
            # Rounded energy values must reproduce the Ritz wavenumber to their quoted precision.
            assert abs(1e8 / wl - (float(r["Ek(cm-1)"]) - float(r["Ei(cm-1)"]))) < 0.01
            selected.append({
                "ion": ion, "label": f"{ion} {target}",
                "ritz_vacuum_wavelength_A": r["ritz_wl_vac(A)"],
                "ritz_wavelength_uncertainty_1sigma_A": err,
                "rest_wavelength_uncertainty_equivalent_m_s":
                    f"{299792458.0 * float(err) / wl:.6f}" if err else "",
                "observed_lab_vacuum_wavelength_A": r.get("obs_wl_vac(A)", ""),
                "observed_lab_wavelength_uncertainty_1sigma_A": r.get("unc_obs_wl", ""),
                "lower_energy_cm_inverse": r["Ei(cm-1)"],
                "upper_energy_cm_inverse": r["Ek(cm-1)"],
                "lower_configuration": r["conf_i"], "lower_term": r["term_i"], "lower_J": r["J_i"],
                "upper_configuration": r["conf_k"], "upper_term": r["term_k"], "upper_J": r["J_k"],
                "transition_type": "E1", "line_reference_code": r.get("line_ref", ""),
                "notes": "Odd parity marked *; data are line identification/calibration input, not a cosmic measurement"
                if ion == "Fe II" else "No explicit ASD wavelength uncertainty; NOT a ppm-grade rest-frequency reference",
            })
    with (DEST / "selected_transitions.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(selected[0]))
        writer.writeheader()
        writer.writerows(selected)
    manifest = {
        "database": "NIST Atomic Spectra Database",
        "version_from_archived_query": version_match.group(1) if version_match else "not recovered",
        "citation": "Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2024), NIST ASD ver. 5.12",
        "doi": "https://doi.org/10.18434/T4W30F",
        "usage": "Small factual data extract. Cite NIST ASD and original measurement sources. No blanket reuse licence inferred.",
        "units": {"wavelength": "vacuum angstrom, explicitly selected show_av=3", "energy": "cm^-1 (E/hc)",
                  "uncertainty": "one standard deviation when supplied by ASD; missing stays missing"},
        "raw_row_counts": {"Fe II": len(raw_fe), "Mg II": len(raw_mg)},
        "selected_row_count": len(selected),
        "selection": "Allowed E1 transitions from Ei=0; Mg II M2 duplicate row retained in raw source and excluded from selected table",
        "caveats": ["Laboratory wavelength uncertainties may share calibration/energy-level covariance",
                    "Ritz and direct observed laboratory wavelengths are distinct, not independent measurements",
                    "No cosmic response coefficient or alpha sensitivity is encoded in this table"],
        "sources": entries,
    }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"raw_row_counts": manifest["raw_row_counts"], "selected": len(selected),
                      "version": manifest["version_from_archived_query"]}))


if __name__ == "__main__":
    main()
