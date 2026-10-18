#!/usr/bin/env python3
"""Mesh an audited prismatic STEP through its original end section.

CAD preprocessor only. No BREP healing, feature removal, structural solution or
inferred slider joint. A conforming sorted-vertex prism split is used instead
of unconstrained 3D Delaunay, avoiding accidental spanwise tetrahedral slivers.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mesh_step_volume import gmsh, validate_linear_mesh


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('step', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--size', type=float, default=.004)
    p.add_argument('--min-size', type=float, default=.0001)
    p.add_argument('--layers', type=int, default=12)
    p.add_argument('--volume-audit', type=Path, required=True)
    a = p.parse_args()
    if not 0 < a.min_size <= a.size or a.layers < 2:
        p.error('Positive mesh sizes and at least two span layers required.')
    source = a.step.resolve(strict=True); source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    audit = json.loads(a.volume_audit.read_text())
    if audit['source_sha256'] != source_hash or not audit['brep_valid'] or audit['units'] != 'm':
        raise ValueError('STEP and audited BREP do not match.')
    gmsh.initialize(readConfigFiles=False)
    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.option.setNumber('General.NumThreads', 1)
        gmsh.option.setString('Geometry.OCCTargetUnit', 'M')
        gmsh.model.occ.importShapes(str(source)); gmsh.model.occ.synchronize()
        volumes = gmsh.model.getEntities(3)
        if len(volumes) != 1:
            raise ValueError('This preprocessor requires one prismatic solid.')
        faces = [tag for _, tag in gmsh.model.getEntities(2)]
        boxes = {tag: np.array(gmsh.model.getBoundingBox(2, tag)) for tag in faces}
        z0 = min(b[2] for b in boxes.values()); z1 = max(b[5] for b in boxes.values())
        caps = [tag for tag, b in boxes.items() if b[5]-b[2] < 1e-5]
        if len(caps) != 2:
            raise ValueError('Exactly two planar span-end faces required.')
        bottom, top = sorted(caps, key=lambda tag: boxes[tag][2])
        if any(gmsh.model.getType(2, tag) != 'Plane' for tag in caps):
            raise ValueError('Span ends must be planes.')
        # Exact cap locations are obtained from their planar CAD parameterization.
        span_values = []
        for tag in (bottom, top):
            lo, hi = gmsh.model.getParametrizationBounds(2, tag)
            span_values.append(gmsh.model.getValue(2, tag, ((np.array(lo)+hi)/2).tolist())[2])
        z0, z1 = span_values; span = z1-z0
        max_normal_span = 0.; samples = 0
        for tag in set(faces)-set(caps):
            lo, hi = map(np.array, gmsh.model.getParametrizationBounds(2, tag))
            uv = np.array([lo+(hi-lo)*[u,v] for u in (.2,.5,.8) for v in (.2,.5,.8)])
            n = np.array(gmsh.model.getNormal(tag, uv.ravel())).reshape(-1,3)
            max_normal_span = max(max_normal_span, float(abs(n[:,2]).max())); samples += len(uv)
            if boxes[tag][2] > z0+1e-5 or boxes[tag][5] < z1-1e-5:
                raise ValueError(f'Face {tag} does not span the complete extrusion.')
        if max_normal_span > 1e-8:
            raise ValueError('CAD is not span-invariant at the geometric control points.')
        gmsh.option.setNumber('Mesh.MeshSizeMax', a.size)
        gmsh.option.setNumber('Mesh.MeshSizeMin', a.min_size)
        gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 24)
        gmsh.option.setNumber('Mesh.RandomSeed', 1)
        gmsh.option.setNumber('Mesh.Algorithm', 6)
        gmsh.model.mesh.generate(2)
        tags, coords, _ = gmsh.model.mesh.getNodes(); original = np.array(coords).reshape(-1,3)
        index = {int(t): i for i,t in enumerate(tags)}
        types, _, blocks = gmsh.model.mesh.getElements(2, bottom)
        tri = np.array([[index[int(t)] for t in row] for typ, block in zip(types,blocks)
                        if typ == 2 for row in np.array(block).reshape(-1,3)], dtype=int)
        used = np.unique(tri); lookup = np.full(len(original),-1,dtype=int); lookup[used] = np.arange(len(used))
        x = original[used]; tri = lookup[tri]; n = len(x)
        if np.max(abs(x[:,2]-z0)) > 1e-9:
            raise ValueError('Cap nodes do not lie on the expected span-end plane.')
        lines = []; line_faces = []
        for _, curve in gmsh.model.getBoundary([(2,bottom)], oriented=False):
            adjacent = [int(tag) for tag in gmsh.model.getAdjacencies(1,curve)[0] if int(tag) != bottom]
            if len(adjacent) != 1 or adjacent[0] in caps:
                raise ValueError('Ambiguous side-face ownership at an end-section curve.')
            types,_,blocks = gmsh.model.mesh.getElements(1,curve)
            for typ,block in zip(types,blocks):
                if typ != 1: raise ValueError('Only straight line segments expected.')
                for row in np.array(block).reshape(-1,2):
                    ids = [lookup[index[int(t)]] for t in row]
                    if min(ids) < 0: raise ValueError('Cap boundary/mesh connectivity mismatch.')
                    lines.append(ids); line_faces.append(adjacent[0])
        lines = np.sort(np.array(lines),axis=1); line_faces = np.array(line_faces)
        nodes = np.vstack([x+[0,0,z-z0] for z in np.linspace(z0,z1,a.layers+1)])
        base = np.sort(tri,axis=1); tetra = []
        for layer in range(a.layers):
            aa,bb,cc = (base+layer*n).T; A,B,C = aa+n,bb+n,cc+n
            tetra.extend([np.stack([aa,bb,cc,C],axis=1),np.stack([aa,bb,B,C],axis=1),np.stack([aa,A,B,C],axis=1)])
        tetra = np.vstack(tetra)
        surface = [tri,tri+a.layers*n]; labels = [np.full(len(tri),bottom),np.full(len(tri),top)]
        for layer in range(a.layers):
            aa,bb = (lines+layer*n).T; A,B = aa+n,bb+n
            surface += [np.stack([aa,bb,B],axis=1),np.stack([aa,B,A],axis=1)]
            labels += [line_faces,line_faces]
        surface = np.vstack(surface); labels = np.concatenate(labels)
        xx = nodes[tetra]; J = np.stack([xx[:,1]-xx[:,0],xx[:,2]-xx[:,0],xx[:,3]-xx[:,0]],axis=2)
        det = np.linalg.det(J); negative = det < 0
        tetra[negative,1],tetra[negative,2] = tetra[negative,2].copy(),tetra[negative,1].copy()
        min_det = validate_linear_mesh(nodes,tetra,surface)
        singular = np.linalg.svd(J,compute_uv=False); ratio = singular[:,-1]/singular[:,0]
        volume = float(abs(det).sum()/6); error = abs(volume/audit['cad_volume_m3']-1)
        if error >= .01: raise ValueError('Extruded mesh exceeds audited volume tolerance.')
        out = a.output; out.mkdir(parents=True,exist_ok=True)
        np.savetxt(out/'nodes.csv',nodes,delimiter=',',header='x_m,y_m,z_m',comments='')
        np.savetxt(out/'tetra.csv',np.column_stack([tetra+1,np.ones(len(tetra),int)]),delimiter=',',fmt='%d',header='n1,n2,n3,n4,cad_volume',comments='')
        np.savetxt(out/'faces.csv',np.column_stack([surface+1,labels]),delimiter=',',fmt='%d',header='n1,n2,n3,cad_face',comments='')
        info = dict(source_name=source.name,source_sha256=source_hash,units='m',gmsh_version=gmsh.__version__,node_count=len(nodes),element_count=len(tetra),face_count=len(surface),cad_volume_count=1,target_size_m=a.size,minimum_size_m=a.min_size,span_layers=a.layers,total_volume_m3=volume,cad_volume_m3=audit['cad_volume_m3'],relative_volume_error=error,min_jacobian_singular_ratio=float(ratio.min()),min_normalized_determinant=min_det,elements_jacobian_ratio_below_1e_4=int((ratio<1e-4).sum()),topology_and_boundary_validated=True,method='original planar CAD end section; conforming 3-tetrahedra prism split; no topology simplification',extrusion_evidence=dict(bottom_face=bottom,top_face=top,span_m=span,side_normal_samples=samples,maximum_span_normal_component=max_normal_span),section_nodes=len(x),section_triangles=len(tri),interpretation='single monolithic CAD solid; no inferred slider joints')
        (out/'mesh_info.json').write_text(json.dumps(info,indent=2)+'\n'); print(json.dumps(info,indent=2),flush=True)
    finally:
        gmsh.finalize()

if __name__ == '__main__': main()
