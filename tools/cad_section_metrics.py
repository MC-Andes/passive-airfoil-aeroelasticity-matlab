#!/usr/bin/env python3
"""Normal material thickness of the unchanged prismatic CAD mesh section."""
import argparse,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser(description=__doc__);p.add_argument('mesh',type=Path);p.add_argument('case',type=Path);p.add_argument('output',type=Path);a=p.parse_args();c=json.loads(a.case.read_text());X=np.loadtxt(a.mesh/'nodes.csv',delimiter=',',skiprows=1);F=np.loadtxt(a.mesh/'faces.csv',delimiter=',',skiprows=1,dtype=int);z0=X[:,2].min();Q=X[:,[0,1]]-np.array([c['geometry']['origin_m'][0],0]);segments=[];labels=[];endpoints=[]
cap=F[F[:,-1]==587,:3]-1;tri=Q[cap];lo=tri.min(axis=1);hi=tri.max(axis=1)
def inside(point):
 use=np.flatnonzero(np.all(point>=lo-1e-13,axis=1)&np.all(point<=hi+1e-13,axis=1));v=tri[use,1]-tri[use,0];w=tri[use,2]-tri[use,0];r=point-tri[use,0];det=v[:,0]*w[:,1]-v[:,1]*w[:,0];b=(r[:,0]*w[:,1]-r[:,1]*w[:,0])/det;d=(v[:,0]*r[:,1]-v[:,1]*r[:,0])/det;return np.any((b>=-1e-10)&(d>=-1e-10)&(b+d<=1+1e-10))
for f in F:
 if f[-1] in (586,587):continue
 ids=f[:3]-1;edge=ids[abs(X[ids,2]-z0)<1e-10]
 if len(edge)==2:segments.append(Q[edge]);labels.append(f[-1]);endpoints.append(edge.tolist())
segments=np.array(segments);labels=np.array(labels);A=segments[:,0];T=segments[:,1]-A
cross=lambda u,v:u[...,0]*v[...,1]-u[...,1]*v[...,0]
rows=[]
for k,(seg,label) in enumerate(zip(segments,labels)):
 point=seg.mean(axis=0);length=np.linalg.norm(T[k]);normal=np.array([-T[k,1],T[k,0]])/length
 if not inside(point+1e-8*normal):normal=-normal
 assert inside(point+1e-8*normal),f'No material side for boundary edge {k}'
 det=cross(normal,T);r=A-point
 with np.errstate(divide='ignore',invalid='ignore'):distance=cross(r,T)/det;fraction=cross(r,normal)/det
 good=(abs(det)>1e-16)&(distance>1e-8)&(fraction>=-1e-9)&(fraction<=1+1e-9)
 assert good.any();hit=np.argmin(np.where(good,distance,np.inf));group=0
 if label in c['regions']['upperTailCadFaces']:group=1
 elif label in c['regions']['lowerTailCadFaces']:group=2
 elif label in [7,9]:group=3
 rows.append([int(label),*point,*normal,float(distance[hit]),int(labels[hit]),length,group])
out=a.output;out.mkdir(parents=True,exist_ok=True);np.savetxt(out/'normal_thickness.csv',rows,delimiter=',',header='cad_face,x_m,z_m,inward_nx,inward_nz,normal_material_thickness_m,opposite_cad_face,segment_length_m,group_upper1_lower2_bar3',comments='');rows=np.array(rows);summary={}
for name,group in [('upper_tail',1),('lower_tail',2),('candidate_bar',3)]:
 r=rows[rows[:,-1]==group];summary[name]=dict(minimum_m=float(r[:,5].min()),median_m=float(np.median(r[:,5])),maximum_m=float(r[:,5].max()),segments=len(r))
summary.update(method='normal ray to first material exit on the original cap mesh; no outer-envelope shortcut',source_sha256=c['geometry']['sourceSha256'],limitations='faceted section; corner normals discontinuous; thick ribs are not assigned a skin thickness')
(out/'normal_thickness.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

# Connected boundary loops identify actual open spanwise cavity ports.
graph={}
for edge,(u,v) in enumerate(endpoints):
 graph.setdefault(u,[]).append(edge);graph.setdefault(v,[]).append(edge)
assert all(len(edges)==2 for edges in graph.values()),'Section boundary must consist of closed nonbranching loops'
unused=set(range(len(endpoints)));loops=[]
while unused:
 edge=min(unused);start=endpoints[edge][0];current=start;nodes=[];owned=[]
 while True:
  nodes.append(current);owned.append(edge);unused.remove(edge);u,v=endpoints[edge];current=v if current==u else u
  if current==start:break
  edge=next(e for e in graph[current] if e!=edge)
  assert edge in unused,'Boundary loop intersects an already visited loop'
 xy=Q[nodes];area=.5*abs(np.sum(cross(xy,np.roll(xy,-1,axis=0))));perimeter=sum(np.linalg.norm(T[e]) for e in owned)
 loops.append(dict(loop=len(loops)+1,area_m2=float(area),perimeter_m=float(perimeter),hydraulic_diameter_m=float(4*area/perimeter),cad_faces=sorted(set(int(labels[e]) for e in owned))))
outer=max(range(len(loops)),key=lambda i:loops[i]['area_m2'])
for i,loop in enumerate(loops):loop['kind']='outer' if i==outer else 'open_spanwise_cavity'
report=dict(source_sha256=c['geometry']['sourceSha256'],span_m=c['geometry']['span_m'],loops=loops,interpretation='ports inferred from connected section boundary loops of the unchanged prismatic solid; actual vent impedance/fixture not measured')
(out/'cavity_ports.json').write_text(json.dumps(report,indent=2)+'\n')
print('Cavity ports:',len(loops)-1,'open area per span end:',sum(x['area_m2'] for x in loops if x['kind']!='outer'))
