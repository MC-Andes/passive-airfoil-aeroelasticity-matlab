#!/usr/bin/env python3
"""Independently audit a STEP BREP using adaptive OpenCASCADE integration.

Optional preprocessing dependency: cadquery-ocp. No structural solver is used.
OpenCASCADE's nonadaptive integration can be inaccurate for trimmed splines.
"""
import argparse
import hashlib
import importlib.metadata
import json
import sys
from pathlib import Path

local = Path(__file__).resolve().parents[1] / ".cache" / "ocp-python"
abi = f"cpython-{sys.version_info.major}{sys.version_info.minor}"
if local.exists() and any((local / "OCP").glob(f"*{abi}*")):
    sys.path.insert(0, str(local))
from OCP.STEPControl import STEPControl_Reader
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.IFSelect import IFSelect_RetDone


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    reader = STEPControl_Reader()
    if reader.ReadFile(str(args.step)) != IFSelect_RetDone:
        raise ValueError("STEP reading failed.")
    reader.SetSystemLengthUnit(1.0)  # OpenCASCADE unit millimetre.
    reader.TransferRoots()
    shape = reader.OneShape()
    valid = BRepCheck_Analyzer(shape).IsValid()
    quick, adaptive = GProp_GProps(), GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, quick)
    error = BRepGProp.VolumeProperties_s(shape, adaptive, 1e-8, True)
    data = dict(source_name=args.step.name,
                source_sha256=hashlib.sha256(args.step.read_bytes()).hexdigest(),
                brep_valid=valid, units="m", cad_volume_m3=adaptive.Mass()/1e9,
                nonadaptive_volume_m3=quick.Mass()/1e9,
                adaptive_tolerance=1e-8, adaptive_reported_relative_error=error,
                method="OpenCASCADE BRepGProp adaptive Gauss integration, OnlyClosed=True",
                wrapper="cadquery-ocp "+importlib.metadata.version("cadquery-ocp"))
    if not valid or data["cad_volume_m3"] <= 0:
        raise ValueError("Invalid or nonpositive CAD solid.")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2)+"\n")
    print(json.dumps(data, indent=2))


if __name__ == "__main__":
    main()
