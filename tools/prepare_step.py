#!/usr/bin/env python3
"""Tessellate STEP in metres and sample outer sections; no material inference.

Requires gmsh and numpy. The OpenCASCADE kernel reads the CAD; Cartesian
control points in STEP are NOT used as an approximation of its boundary.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

LOCAL_PACKAGES = Path(__file__).resolve().parents[1] / ".cache" / "cad-python"
if LOCAL_PACKAGES.exists():
    sys.path.insert(0, str(LOCAL_PACKAGES))
import gmsh


def outer_section(nodes, triangles, axis, location, chord_axis, height_axis, count):
    """Intersect boundary triangles with a plane and retain outer envelopes.

    Appropriate only for sections representable by single-valued upper/lower
    curves. Internal holes are deliberately omitted from the aerodynamic skin.
    """
    segments = []
    scale = np.ptp(nodes[:, axis])
    tol = max(scale, 1.0) * 1e-10
    for tri in nodes[triangles]:
        d = tri[:, axis] - location
        if np.all(d > tol) or np.all(d < -tol) or np.all(abs(d) < tol):
            continue
        intersections = []
        for i, j in ((0, 1), (1, 2), (2, 0)):
            if abs(d[i]) <= tol:
                intersections.append(tri[i])
            if d[i] * d[j] < 0:
                intersections.append(tri[i] + (tri[j] - tri[i]) * d[i] / (d[i] - d[j]))
        if len(intersections) < 2:
            continue
        pts = np.unique(np.round(intersections, 12), axis=0)
        if len(pts) == 2 and np.linalg.norm(pts[1] - pts[0]) > tol:
            segments.append(pts[:, [chord_axis, height_axis]])
    if not segments:
        raise ValueError(f"No section at coordinate {location:g}.")
    segments = np.asarray(segments)
    xmin, xmax = segments[:, :, 0].min(), segments[:, :, 0].max()
    xi = (1 - np.cos(np.linspace(0, np.pi, count))) / 2
    x = xmin + (xmax - xmin) * xi
    upper, lower = [], []
    for xx in x:
        heights = []
        for a, b in segments:
            if min(a[0], b[0]) - tol <= xx <= max(a[0], b[0]) + tol:
                dx = b[0] - a[0]
                if abs(dx) <= tol:
                    heights.extend([a[1], b[1]])
                else:
                    ratio = np.clip((xx - a[0]) / dx, 0, 1)
                    heights.append(a[1] + ratio * (b[1] - a[1]))
        if not heights:
            raise ValueError("Section has gaps: automatic outer-envelope extraction is unsuitable.")
        lower.append(min(heights))
        upper.append(max(heights))
    lower, upper = np.asarray(lower), np.asarray(upper)
    return np.column_stack([xi, x, upper, lower, (upper + lower) / 2]), segments


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--axes", default="0,2,1", help="CAD chord,span,height axes (0,1,2)")
    parser.add_argument("--mesh-size", type=float, default=0.006, help="Surface edge target in metres")
    parser.add_argument("--sections", type=int, default=5)
    parser.add_argument("--samples", type=int, default=161)
    args = parser.parse_args()
    axes = list(map(int, args.axes.split(",")))
    if sorted(axes) != [0, 1, 2] or args.mesh_size <= 0 or args.sections < 2 or args.samples < 3:
        parser.error("Invalid axes, mesh size, or sample count.")
    step = args.step.resolve(strict=True)
    args.output.mkdir(parents=True, exist_ok=True)
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setString("Geometry.OCCTargetUnit", "M")
        gmsh.model.occ.importShapes(str(step))
        gmsh.model.occ.synchronize()
        solids = gmsh.model.getEntities(3)
        if not solids:
            raise ValueError("No solid BREP found in STEP.")
        bounds = np.array([gmsh.model.getBoundingBox(*item) for item in solids])
        kernel_minimum, kernel_maximum = bounds[:, :3].min(axis=0), bounds[:, 3:].max(axis=0)
        gmsh.option.setNumber("Mesh.MeshSizeMax", args.mesh_size)
        gmsh.option.setNumber("Mesh.MeshSizeMin", args.mesh_size / 10)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 12)
        gmsh.model.mesh.generate(2)
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        nodes = np.asarray(xyz).reshape(-1, 3)
        minimum, maximum = nodes.min(axis=0), nodes.max(axis=0)
        tag_index = {int(tag): i for i, tag in enumerate(tags)}
        types, _, connectivity = gmsh.model.mesh.getElements(2)
        triangle_blocks = []
        for kind, block in zip(types, connectivity):
            name, _, _, number_nodes, _, _ = gmsh.model.mesh.getElementProperties(kind)
            if name.startswith("Triangle"):
                block = np.asarray(block).reshape(-1, number_nodes)[:, :3]
                triangle_blocks.append(np.array([[tag_index[int(t)] for t in row] for row in block]))
        if not triangle_blocks:
            raise ValueError("No triangular surface mesh generated.")
        triangles = np.vstack(triangle_blocks)
        np.savetxt(args.output / "cad_nodes.csv", nodes, delimiter=",", header="x_m,y_m,z_m", comments="")
        np.savetxt(args.output / "cad_triangles.csv", triangles + 1, delimiter=",", fmt="%d", header="i,j,k", comments="")
        gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.write(str(args.output / "cad_surface.stl"))
        chord_axis, span_axis, height_axis = axes
        # Interior planes avoid coplanar end-cap triangles. Edge stations will
        # be extrapolated only if a user explicitly requests a surface model.
        fractions = np.linspace(0.01, 0.99, args.sections)
        records = []
        for i, fraction in enumerate(fractions):
            location = minimum[span_axis] + fraction * (maximum[span_axis] - minimum[span_axis])
            section, segments = outer_section(nodes, triangles, span_axis, location, chord_axis, height_axis, args.samples)
            filename = f"section_{i + 1:02d}.csv"
            np.savetxt(args.output / filename, section, delimiter=",", header="xi,x_m,upper_m,lower_m,camber_m", comments="")
            np.savetxt(args.output / f"section_segments_{i + 1:02d}.csv", segments.reshape(-1, 4), delimiter=",", header="x1_m,z1_m,x2_m,z2_m", comments="")
            records.append({"fraction": float(fraction), "span_m": float(location), "file": filename,
                            "chord_m": float(np.ptp(section[:, 1]))})
        relative_source = None
        try:
            relative_source = str(step.relative_to(args.output.resolve().parents[1]))
        except ValueError:
            pass
        info = {"source_name": step.name, "source_relative_path": relative_source,
                "source_sha256": hashlib.sha256(step.read_bytes()).hexdigest(), "units": "m",
                "kernel": "OpenCASCADE via Gmsh", "gmsh_version": gmsh.__version__,
                "cad_axes_zero_based": axes, "solid_count": len(solids),
                "face_count": len(gmsh.model.getEntities(2)), "node_count": len(nodes),
                "triangle_count": len(triangles), "minimum_m": minimum.tolist(),
                "maximum_m": maximum.tolist(), "sections": records,
                "bounding_box_source": "tessellated boundary nodes",
                "kernel_minimum_m": kernel_minimum.tolist(), "kernel_maximum_m": kernel_maximum.tolist(),
                "interpretation": "Outer envelopes only; no materials, constraints, joints or stiffness inferred. Inspect before use."}
        (args.output / "cad_info.json").write_text(json.dumps(info, indent=2) + "\n")
        print(json.dumps({key: info[key] for key in ("units", "solid_count", "face_count", "node_count", "triangle_count", "minimum_m", "maximum_m")}, indent=2))
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    main()
