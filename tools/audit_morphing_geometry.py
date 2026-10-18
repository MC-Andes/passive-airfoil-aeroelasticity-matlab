#!/usr/bin/env python3
"""Verify source-bound socket/slider signatures on the updated STEP.

Preprocessing only: uses Gmsh/OpenCASCADE, never a structural/CFD solver.
Run audit_step_volume.py first. Inventor IPT is recorded by hash, not solved.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from mesh_step_volume import gmsh


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('step',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--case',type=Path,default=Path(__file__).resolve().parents[1]/'cases/morphing_wing/case.json')
    p.add_argument('--ipt',type=Path)
    a=p.parse_args();c=json.loads(a.case.read_text());root=Path(__file__).resolve().parents[1]
    reference=json.loads((root/c['geometry']['audit']).read_text());digest=hashlib.sha256(a.step.read_bytes()).hexdigest()
    if digest!=c['geometry']['sourceSha256']:raise ValueError('Unbound STEP: identify interfaces again; do not reuse CAD labels.')
    gmsh.initialize(readConfigFiles=False)
    try:
        gmsh.option.setNumber('General.Terminal',0);gmsh.option.setString('Geometry.OCCTargetUnit','M')
        gmsh.model.occ.importShapes(str(a.step.resolve()));gmsh.model.occ.synchronize()
        faces=[]
        for _,tag in gmsh.model.getEntities(2):
            lo,hi=map(np.asarray,gmsh.model.getParametrizationBounds(2,tag))
            faces.append(dict(tag=tag,type=gmsh.model.getType(2,tag),box=gmsh.model.getBoundingBox(2,tag),center=gmsh.model.occ.getCenterOfMass(2,tag),normal=gmsh.model.getNormal(tag,((lo+hi)/2).tolist()).tolist(),area_m2=gmsh.model.occ.getMass(2,tag)))
        by={f['tag']:f for f in faces};checks=[]
        for group in ('supports','slider'):
            for tag in c[group]['cadFaces']:
                expected=reference['selected_faces'][str(tag)];found=by[tag]
                assert found['type']==expected['type'] and np.linalg.norm(np.array(found['center'])-expected['center'])<1e-10
                assert abs(found['area_m2']/expected['area_m2']-1)<1e-9
                checks.append(tag)
        span=float(max(f['box'][5] for f in faces)-min(f['box'][2] for f in faces));partial=[f['tag'] for f in faces if f['box'][5]-f['box'][2]>.001 and f['box'][5]-f['box'][2]<.35]
        report=dict(source_sha256=digest,solids=len(gmsh.model.getEntities(3)),faces=len(faces),curves=len(gmsh.model.getEntities(1)),verified_interface_faces=checks,nonprismatic_partial_span_faces=partial,full_inventory=faces,approximate_kernel_span_m=span,scope='signature audit only; origin and nominal dimensions from declared study; no automatic discovery of unconfirmed contact or printing law')
        if a.ipt:report['ipt_sha256']=hashlib.sha256(a.ipt.read_bytes()).hexdigest();report['ipt_equivalence']='native file preserved; Inventor feature-tree equivalence not verified by this STEP audit'
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2)+'\n')
        print(f'PASS: {report["solids"]} solid, {len(faces)} faces, interfaces {checks}; nonprismatic CAD; source {digest}')
    finally:gmsh.finalize()


if __name__=='__main__':main()
