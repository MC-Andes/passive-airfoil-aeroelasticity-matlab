#!/usr/bin/env python3
"""CAD preprocessing only: STEP -> linear tetrahedra for the MATLAB solver.

Gmsh/OpenCASCADE generate the mesh, never stiffness, loads or a solution.
Coordinates remain in original CAD axes, in metres. CAD faces/volumes are
exported as labels; these labels do not identify materials, joints or supports.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

local = Path(__file__).resolve().parents[1] / ".cache" / "cad-python"
if local.exists():
    sys.path.insert(0, str(local))
import gmsh



def validate_linear_mesh(nodes, tetra, triangles):
    """Reject invalid exports; never remove CAD features/elements to force a solve."""
    corners = nodes[tetra]
    jac = np.stack([corners[:, 1]-corners[:, 0], corners[:, 2]-corners[:, 0],
                    corners[:, 3]-corners[:, 0]], axis=2)
    scale = np.linalg.norm(jac, axis=1).max(axis=1)
    normalized = abs(np.linalg.det(jac))/scale**3
    if np.any(normalized <= 1e-13):
        raise ValueError(f"Degenerate volume mesh: {int((normalized <= 1e-13).sum())} "
                         "elements fail the MATLAB determinant threshold; remesh without deleting features.")
    faces = np.sort(np.concatenate([tetra[:, face] for face in
                                   ([1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1])]), axis=1)
    unique, counts = np.unique(faces, axis=0, return_counts=True)
    if counts.max() > 2:
        raise ValueError("Nonmanifold/duplicate tetrahedral faces; mesh export rejected.")
    surface = np.sort(triangles, axis=1)
    boundary = unique[counts == 1]
    if (len(surface) != len(boundary) or len(np.unique(surface, axis=0)) != len(surface)
            or len(np.unique(np.concatenate([surface, boundary]), axis=0)) != len(surface)):
        raise ValueError("CAD-labeled triangles do not match the volume boundary; mesh export rejected.")
    return float(normalized.min())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--size", type=float, default=0.006)
    parser.add_argument("--min-size", type=float, default=0.0002)
    parser.add_argument("--verbose", action="store_true", help="Print Gmsh diagnostics")
    parser.add_argument("--threads", type=int, default=1, help="Gmsh worker count (default 1)")
    parser.add_argument("--netgen", action="store_true", help="Optional Netgen optimization")
    parser.add_argument("--volume-audit", type=Path, help="JSON from audit_step_volume.py")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("threads must be positive")
    if not 0 < args.min_size <= args.size:
        parser.error("Require 0 < min-size <= size, in metres.")
    source = args.step.resolve(strict=True)
    args.output.mkdir(parents=True, exist_ok=True)
    gmsh.initialize(readConfigFiles=False)
    try:
        gmsh.option.setNumber("General.Terminal", int(args.verbose))
        gmsh.option.setNumber("General.NumThreads", args.threads)
        gmsh.option.setString("Geometry.OCCTargetUnit", "M")
        gmsh.model.occ.importShapes(str(source))
        gmsh.model.occ.synchronize()
        solids = gmsh.model.getEntities(3)
        if not solids:
            raise ValueError("STEP has no solid volumes.")
        gmsh.option.setNumber("Mesh.MeshSizeMax", args.size)
        gmsh.option.setNumber("Mesh.MeshSizeMin", args.min_size)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 16)
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.RandomSeed", 1)
        print("Generating STEP volume mesh in metres...", flush=True)
        gmsh.model.mesh.generate(3)
        if args.netgen:
            gmsh.model.mesh.optimize("Netgen")
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        original = np.asarray(xyz).reshape(-1, 3)
        index = {int(t): i for i, t in enumerate(tags)}
        tetra, volume_id = [], []
        for _, tag in solids:
            for kind, block in zip(gmsh.model.mesh.getElements(3, tag)[0],
                                   gmsh.model.mesh.getElements(3, tag)[2]):
                if kind != 4:
                    raise ValueError("Only linear four-node tetrahedra supported.")
                block = np.asarray(block).reshape(-1, 4)
                tetra.extend([[index[int(t)] for t in row] for row in block])
                volume_id.extend([tag] * len(block))
        tetra = np.asarray(tetra, dtype=int)
        used = np.unique(tetra)
        lookup = np.full(len(original), -1, dtype=int)
        lookup[used] = np.arange(len(used))
        nodes, tetra = original[used], lookup[tetra]
        triangles, face_id = [], []
        for _, tag in gmsh.model.getEntities(2):
            types, _, blocks = gmsh.model.mesh.getElements(2, tag)
            for kind, block in zip(types, blocks):
                if kind != 2:
                    raise ValueError("Only linear triangles supported.")
                rows = [[lookup[index[int(t)]] for t in row]
                        for row in np.asarray(block).reshape(-1, 3)]
                if np.any(np.asarray(rows) < 0):
                    raise ValueError("Surface nodes missing from volume mesh.")
                triangles.extend(rows)
                face_id.extend([tag] * len(rows))
        triangles = np.asarray(triangles, dtype=int)
        # Orient tetrahedra positively in the original coordinate system.
        xyz = nodes[tetra]
        det = np.linalg.det(np.stack([xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0],
                                     xyz[:, 3]-xyz[:, 0]], axis=2))
        negative = det < 0
        tetra[negative, 1], tetra[negative, 2] = (tetra[negative, 2].copy(),
                                                 tetra[negative, 1].copy())
        min_normalized_det = validate_linear_mesh(nodes, tetra, triangles)
        np.savetxt(args.output / "nodes.csv", nodes, delimiter=",",
                   header="x_m,y_m,z_m", comments="")
        np.savetxt(args.output / "tetra.csv", np.column_stack([tetra+1, volume_id]),
                   delimiter=",", fmt="%d", header="n1,n2,n3,n4,cad_volume", comments="")
        np.savetxt(args.output / "faces.csv", np.column_stack([triangles+1, face_id]),
                   delimiter=",", fmt="%d", header="n1,n2,n3,cad_face", comments="")
        quick_volume = sum(gmsh.model.occ.getMass(dim, tag) for dim, tag in solids)
        cad_volume = quick_volume
        volume_method = "Gmsh nonadaptive OCC integration; not independently audited"
        if args.volume_audit:
            audit = json.loads(args.volume_audit.read_text())
            if (audit["source_sha256"] != hashlib.sha256(source.read_bytes()).hexdigest()
                    or audit["units"] != "m" or not audit["brep_valid"]):
                raise ValueError("Volume audit does not match this STEP solid.")
            cad_volume = audit["cad_volume_m3"]
            volume_method = audit["method"]
        total_volume = float(np.sum(abs(det))/6)
        singular = np.linalg.svd(np.stack([xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0],
                                         xyz[:, 3]-xyz[:, 0]], axis=2), compute_uv=False)
        ratio = singular[:, -1]/singular[:, 0]
        info = dict(source_name=source.name,
                    source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                    units="m", gmsh_version=gmsh.__version__,
                    node_count=len(nodes), element_count=len(tetra),
                    face_count=len(triangles), cad_volume_count=len(solids),
                    target_size_m=args.size, minimum_size_m=args.min_size,
                    threads=args.threads, netgen_optimization=args.netgen,
                    algorithm_3d=1, random_seed=1, curvature_sizing=16,
                    user_config_files_read=False,
                    total_volume_m3=total_volume, cad_volume_m3=cad_volume,
                    cad_volume_nonadaptive_m3=quick_volume,
                    volume_reference_method=volume_method,
                    relative_volume_error=abs(total_volume-cad_volume)/cad_volume,
                    min_jacobian_singular_ratio=float(ratio.min()),
                    min_normalized_determinant=min_normalized_det,
                    topology_and_boundary_validated=True,
                    elements_jacobian_ratio_below_1e_4=int((ratio < 1e-4).sum()),
                    status="PRELIMINARY mesh; perform mesh convergence before quantitative prediction",
                    interpretation="One connected CAD solid; no slider kinematics or material inferred.")
        (args.output / "mesh_info.json").write_text(json.dumps(info, indent=2)+"\n")
        print(json.dumps(info, indent=2))
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    main()
