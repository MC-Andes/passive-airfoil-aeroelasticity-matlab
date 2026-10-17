#!/usr/bin/env python3
"""Evaluate CAD-projected quadratic midside nodes without altering the BREP.

Preprocessing experiment. A folded projection is rejected and reported, never
silently blended back. Structural FE evaluation remains in MATLAB.
"""
import argparse,json,hashlib,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from mesh_step_volume import gmsh
p=argparse.ArgumentParser(description=__doc__);p.add_argument('step',type=Path);p.add_argument('mesh',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
info=json.loads((a.mesh/'mesh_info.json').read_text());assert hashlib.sha256(a.step.read_bytes()).hexdigest()==info['source_sha256']
X=np.loadtxt(a.mesh/'nodes.csv',delimiter=',',skiprows=1);T=np.loadtxt(a.mesh/'tetra.csv',delimiter=',',skiprows=1,dtype=int);F=np.loadtxt(a.mesh/'faces.csv',delimiter=',',skiprows=1,dtype=int);T4=T[:,:4]-1;F3=F[:,:3]-1
pairs=np.array([[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]])
E,inv=np.unique(np.sort(T4[:,pairs].reshape(-1,2),axis=1),axis=0,return_inverse=True);T10=np.column_stack([T4,len(X)+inv.reshape(-1,6)])
mid=(X[E[:,0]]+X[E[:,1]])/2;lookup={tuple(e):i for i,e in enumerate(E)};ownership={}
for face,label in zip(F3,F[:,3]):
 for i,j in [(0,1),(1,2),(2,0)]:ownership.setdefault(lookup[tuple(sorted([face[i],face[j]]))],set()).add(int(label))
gmsh.initialize(readConfigFiles=False)
try:
 gmsh.option.setNumber('General.Terminal',0);gmsh.option.setString('Geometry.OCCTargetUnit','M');gmsh.model.occ.importShapes(str(a.step));gmsh.model.occ.synchronize()
 curves={f:set(tag for dim,tag in gmsh.model.getBoundary([(2,f)],oriented=False)) for f in set(F[:,3])}
 projected=mid.copy();ambiguous=0
 for i,labels in ownership.items():
  common=set.intersection(*(curves[f] for f in labels))
  if len(labels)>1 and common:
   candidates=[np.array(gmsh.model.getClosestPoint(1,c,mid[i].tolist())[0]) for c in common];projected[i]=min(candidates,key=lambda x:np.linalg.norm(x-mid[i]))
  elif len(labels)==1:projected[i]=gmsh.model.getClosestPoint(2,next(iter(labels)),mid[i].tolist())[0]
  else:ambiguous+=1
 shift=np.linalg.norm(projected-mid,axis=1);XX=np.vstack([X,projected]);controls=np.vstack([np.eye(4),np.full((1,4),.25),np.array([[.5,.5,0,0],[0,.5,.5,0],[.5,0,.5,0],[.5,0,0,.5],[0,.5,0,.5],[0,0,.5,.5]])]);gradL=np.vstack([-np.ones(3),np.eye(3)]);fold=np.zeros(len(T4),bool);minratio=1.
 local=XX[T10]-XX[T10[:,0],None,:]
 original=np.linalg.det(np.stack([X[T4[:,i]]-X[T4[:,0]] for i in [1,2,3]],axis=2))
 for L in controls:
  G=np.vstack([(4*L-1)[:,None]*gradL,np.array([4*(L[i]*gradL[j]+L[j]*gradL[i]) for i,j in pairs])]);det=np.linalg.det(np.einsum('eia,ib->eab',local,G));fold|=(det*original<=0);minratio=min(minratio,float(np.min(abs(det/original))))
 quadrature_volume=0.;gauss,weights=np.polynomial.legendre.leggauss(4);gauss=(gauss+1)/2;weights=weights/2
 for i,r in enumerate(gauss):
  for j,ss in enumerate(gauss):
   for k,t in enumerate(gauss):
    L=np.array([1-r-(1-r)*ss-(1-r)*(1-ss)*t,r,(1-r)*ss,(1-r)*(1-ss)*t]);G=np.vstack([(4*L-1)[:,None]*gradL,np.array([4*(L[u]*gradL[v]+L[v]*gradL[u]) for u,v in pairs])]);det=np.linalg.det(np.einsum('eia,ib->eab',local,G));fold|=det*original<=0;quadrature_volume+=float(np.sum(abs(det)))*weights[i]*weights[j]*weights[k]*(1-r)**2*(1-ss)
 info['total_volume_m3']=quadrature_volume;info['relative_volume_error']=abs(quadrature_volume/info['cad_volume_m3']-1)
 out=a.output;out.mkdir(parents=True,exist_ok=True);report=dict(source_sha256=info['source_sha256'],maximum_projection_m=float(shift.max()),rms_projection_m=float(np.sqrt(np.mean(shift**2))),projected_boundary_edges=len(ownership),ambiguous_edges=ambiguous,integrated_curved_volume_m3=quadrature_volume,relative_volume_error=info['relative_volume_error'],folded_elements_at_controls_and_quadrature=int(fold.sum()),minimum_abs_jacobian_ratio=minratio,accepted=not fold.any() and ambiguous==0,scope='quadratic CAD projection evaluated; unchanged original BREP; rejected projections do not replace converged affine geometry')
 (out/'projection_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
 if report['accepted']:
  FE=np.array([[lookup[tuple(sorted([f[i],f[j]]))]+len(X) for i,j in [(0,1),(1,2),(2,0)]] for f in F3]);F6=np.column_stack([F3,FE])
  np.savetxt(out/'nodes.csv',XX,delimiter=',',header='x_m,y_m,z_m',comments='');np.savetxt(out/'tetra.csv',np.column_stack([T10+1,T[:,-1]]),delimiter=',',fmt='%d',header='n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,cad_volume',comments='');np.savetxt(out/'faces.csv',np.column_stack([F6+1,F[:,-1]]),delimiter=',',fmt='%d',header='n1,n2,n3,n4,n5,n6,cad_face',comments='');info['node_count']=len(XX);info['curved_projection']=report;(out/'mesh_info.json').write_text(json.dumps(info,indent=2)+'\n')
finally:gmsh.finalize()
