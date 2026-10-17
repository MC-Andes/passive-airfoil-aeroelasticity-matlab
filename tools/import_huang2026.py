#!/usr/bin/env python3
"""Extract a small, attributed modal benchmark from the original v2 workbook.

Usage: python tools/import_huang2026.py path/to/original.xlsx
Uses only Python's standard library; no workbook modifications or recalculation.
The original workbook is obtained from the URL in the provenance JSON.
"""
import csv
import hashlib
import json
import math
import statistics
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SHA = "34545c9f937404864ad34a4d6dd5c03c0f23b20138115a5bf3dc61fea5749408"
URL = ("https://data.mendeley.com/public-files/datasets/8jxw5533py/files/"
       "d8280c9e-63d7-43fb-8d27-a7e511452f9c/file_downloaded")


def extract(path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != SHA:
        raise ValueError("Workbook differs from the documented dataset v2 input")
    ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(path) as archive:
        strings = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        strings = ["".join(si.itertext()) for si in strings]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    cells = {}
    for c in sheet.findall(".//s:c", ns):
        v = c.find("s:v", ns)
        if v is not None:
            cells[c.attrib["r"]] = strings[int(v.text)] if c.get("t") == "s" else float(v.text)
    # Key cases by geometry, since Test4/Test6 labels disagree with article Table 1.
    groups = [(2, 1.8, 0), (16, 2.5, 0), (31, 5, 0), (44, 5, .25),
              (58, 5, .5), (75, 5, .75), (95, 3, .5), (110, 4, .5), (126, 6, .5)]
    rows = []
    for header, ar, flex in groups:
        case_id = f"HZ_AR{ar:g}_f{flex:g}"
        for rank in range(1, 4):
            r = header + 1 + rank
            values = [cells[f"{c}{r}"] for c in ("D", "H", "L")]
            mean, sd = cells[f"N{r}"], cells[f"O{r}"]
            assert math.isclose(sd, statistics.stdev(values), rel_tol=1e-12)
            rows.append([case_id, ar, flex, rank, cells[f"A{r}"], mean, sd,
                         cells[f"B{r}"], cells[f"F{r}"], cells[f"J{r}"],
                         f"Modal Test Result!N{r}:O{r}", cells[f"A{header}"],
                         statistics.mean(values), int(math.isclose(mean, statistics.mean(values), rel_tol=1e-12))])
    out = ROOT / "data/benchmarks"
    out.mkdir(parents=True, exist_ok=True)
    target = out / "huang2026_modal.csv"
    with target.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["case_id", "aspect_ratio", "flexible_chord_fraction", "frequency_rank",
                         "source_mode_label", "omega_nd_mean", "omega_nd_sample_sd",
                         "replicate_1_Hz", "replicate_2_Hz", "replicate_3_Hz",
                         "source_cells", "source_specimen_label", "all_three_Omega_mean",
                         "published_mean_uses_all_three"])
        writer.writerows(rows)
    provenance = {
        "authors": ["Ke Huang", "Jiaying Zhang"], "dataset_doi": "10.17632/8jxw5533py.2",
        "article_doi": "10.1016/j.dib.2026.113015", "dataset_version": 2,
        "dataset_url": "https://data.mendeley.com/datasets/8jxw5533py/2",
        "workbook_url": URL, "workbook_sha256": SHA,
        "workbook_name": "1-1 Modal Test Experiment and Numerical Results.xlsx",
        "license": "CC BY 4.0", "license_url": "https://creativecommons.org/licenses/by/4.0/",
        "changes": "Selected first three frequency-ordered modes in nine PLA geometry groups; "
                   "cached dimensionless means/sample SD and individual Hz transcribed without recalculation; "
                   "geometry-based IDs and separate all-three means added. Sample SD checked against "
                   "cached replicate values; published subset means preserved and flagged.",
        "nominal_geometry": {"chord_m": .05, "rigid_thickness_m": .004, "flexible_thickness_m": .001,
                             "span_m": "aspect_ratio * chord_m", "flexible_region": "aft chord strip"},
        "normalization": "Omega = 2*pi*f*L^2*sqrt(rho*tR / (E*tR^3/(12*(1-nu^2))))",
        "nominal_model_parameters": {"E_Pa": 2.11e9, "nu": .3, "rho_kg_m3": 1180},
        "model_parameter_source": "AR5fef0d5PPModel.bdf MAT1 and PSHELL in the same dataset",
        "support": "y=0 on rigid chord region only for stepped plates; whole root for uniform",
        "support_source": "AR5fef0d5PPModel.bdf SPC1 set 11: x=0..0.025 m, y=0",
        "cautions": ["Workbook Test4/Test6 geometry labels disagree with article Table 1; keyed by AR/f.",
                     "Published means in rows 46:48 use H/L and rows 97:99 use D/H, while SD uses D/H/L. "
                     "Published means are preserved; a separate all-three mean and flag expose the difference.",
                     "Source normalization uses specimen-specific dimensions and materials; raw Hz are not "
                     "direct targets for a single nominal material.",
                     "Sample SD is replicate spread, not a confidence interval or total uncertainty.",
                     "Only first bending mode has an unambiguous identification in this reduced extract; "
                     "other comparisons use frequency order, not measured mode-shape matching."]}
    (out / "huang2026_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Extracted {len(rows)} observations to {target}")


if __name__ == "__main__":
    extract(Path(sys.argv[1]))
