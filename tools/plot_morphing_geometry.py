#!/usr/bin/env python3
"""Plot actual prepared sections of the updated CAD; no invented geometry."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle

p=argparse.ArgumentParser(description=__doc__);p.add_argument('surface',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
root=Path(__file__).resolve().parents[1];c=json.loads((root/'cases/morphing_wing/case.json').read_text());audit=json.loads((root/c['geometry']['audit']).read_text());info=json.loads((a.surface/'cad_info.json').read_text());assert info['source_sha256']==c['geometry']['sourceSha256']
fig,axs=plt.subplots(3,1,figsize=(10,7),layout='constrained',sharex=True,sharey=True)
for ax,k in zip(axs,[1,3,5]):
 seg=np.loadtxt(a.surface/f'section_segments_{k:02d}.csv',delimiter=',',skiprows=1).reshape(-1,2,2)
 seg=(seg-np.array(c['geometry']['origin_m'])[[0,2]])/c['geometry']['chord_m'];ax.add_collection(LineCollection(seg,colors='#2d5a4e',linewidths=.65))
 h=audit['shaft_socket']['side_m']/2;xy=np.array([.0625-h,.0043-h])/.25;ax.add_patch(Rectangle(xy,2*h/.25,2*h/.25,fill=False,color='#8b3a2f',lw=1.3))
 face=audit['selected_faces']['8'];box=np.array(face['box']);xx=(face['center'][0]-c['geometry']['origin_m'][0])/.25;zz=(box[[1,4]]-c['geometry']['origin_m'][2])/.25;ax.plot([xx,xx],zz,color='#245da8',lw=2)
 ax.set_aspect('equal');ax.set(xlim=(-.025,1.01),ylim=(-.055,.095),ylabel='z/c');ax.text(.01,.9,f"η = {info['sections'][k-1]['fraction']:.2f}",transform=ax.transAxes,fontsize=11)
 ax.spines[['top','right']].set_visible(False)
axs[-1].set_xlabel('x/c, referencia nominal c = 250 mm')
a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=180);plt.close(fig);print(a.output)
