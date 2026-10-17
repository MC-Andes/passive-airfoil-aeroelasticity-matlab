#!/usr/bin/env python3
"""Publish figures from completed local CAD campaigns, never synthetic fields.

Optional plotting dependencies: tools/requirements-portal.txt.
First run export_portal_fields(output) in MATLAB to sample the actual motion.
Large original CAD and campaign files remain outside Git.
"""
import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap
import imageio_ffmpeg
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
GREEN, RED, INK, PAPER = '#2d5a4e', '#8b3a2f', '#1a1a1a', '#faf9f7'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 15, 'axes.labelcolor': INK,
                     'text.color': INK, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': PAPER, 'axes.facecolor': PAPER, 'savefig.facecolor': PAPER})


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''): h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=True)
    sources = []

    def source(value):
        path = ROOT / value
        sources.append({'path': value, 'sha256': digest(path)})
        return path

    def data(value):
        with source(value).open(newline='') as stream: return list(csv.DictReader(stream))

    def save(fig, name):
        fig.savefig(out / name, dpi=180, bbox_inches='tight'); plt.close(fig)

    c = json.loads(source('cases/cad_original/extruded-fine.json').read_text())
    mesh = 'results/cad_original/meshes/extruded_4/'
    X = np.loadtxt(source(mesh+'nodes.csv'), delimiter=',', skiprows=1)
    X = X[:, np.array(c['geometry']['axisOrder'])-1] - c['geometry']['origin_m']
    faces = np.loadtxt(source(mesh+'faces.csv'), delimiter=',', skiprows=1, dtype=int)
    tri = faces[faces[:, -1] == 587, :3]-1
    assert len(tri) > 100, 'Original CAD cap is missing'
    chord = c['geometry']['chord_m']; vertices = X[tri][:, :, [0, 2]] / chord
    fig, ax = plt.subplots(figsize=(7.2, 4))
    ax.add_collection(PolyCollection(vertices, facecolors='#e2e8e3', edgecolors='#6e8278', linewidths=.13))
    ax.axvspan(0, .3, facecolor=RED, alpha=.12, hatch='///', edgecolor=RED, label='Apoyo supuesto: x/c ≤ 0,30')
    ax.autoscale(); ax.set_aspect('equal'); ax.set(xlabel='x/c', ylabel='z/c', xlim=(-.03, 1.03))
    save(fig, 'geometry.png')

    rows = data('results/cad_original/campaign_20261017T190836945/summary.csv')
    dofs = [3*int(r['nodes'])/1e6 for r in rows]
    te = [100*1e3*float(r['maxSensor_m_per_Pa']) for r in rows]
    with (out/'structural-convergence.csv').open('w') as stream:
        stream.write('dofs,TE_100Pa_mm\n')
        for r, value in zip(rows,te): stream.write(f"{3*int(r['nodes'])},{value:.12g}\n")
    fig, ax = plt.subplots(figsize=(7.2,4.2)); ax.plot(dofs, te, 'o-', color=GREEN)
    ax.set(xlabel='Grados de libertad [millones]', ylabel='Máximo TE a 100 Pa [mm]'); ax.grid(alpha=.2)
    save(fig, 'convergence.png')

    review = json.loads(source('results/cad_original/final_review/acceptance.json').read_text())
    source('results/verification/summary.json')
    dynamic = review['manifest']['timeHalf128']
    folder = Path(dynamic); folder = folder if folder.is_absolute() else ROOT/folder
    history = folder/'U_10/te_history.csv'
    sources.append({'path':str(history.relative_to(ROOT)), 'sha256':digest(history)})
    values = np.loadtxt(history, delimiter=','); keep = np.arange(0,len(values),16)
    samples = values[keep][:,[0,5]]; np.savetxt(out/'te-history.csv', samples, delimiter=',', header='time_s,TE_midspan_relative_to_trim_m', comments='')
    fig, ax = plt.subplots(figsize=(7.2,4.2)); ax.plot(samples[:,0],1e6*samples[:,1],color=GREEN,lw=1.5)
    ax.axhline(0,color=INK,lw=.5); ax.set(xlabel='Tiempo físico [s]',ylabel='TE respecto al equilibrio [µm]');ax.grid(alpha=.2)
    save(fig, 'free-response.png')

    rows = data('results/tail_benchmarks/modal_comparison.csv'); rows = [r for r in rows if int(r['frequency_rank'])==1]
    fig, ax = plt.subplots(figsize=(7.2,5.5)); y=np.arange(len(rows))
    ax.errorbar([float(r['experiment_Omega']) for r in rows],y,xerr=[float(r['experiment_sample_sd']) for r in rows],fmt='o',color=INK,label='Ensayo ± SD',capsize=3)
    ax.scatter([float(r['model_Omega']) for r in rows],y,marker='x',s=55,color=RED,label='MITC4')
    ax.set_yticks(y,[r['case_id'].replace('HZ_','') for r in rows],fontsize=13);ax.invert_yaxis()
    ax.set(xlabel='Frecuencia adimensional Ω');ax.grid(axis='x',alpha=.2)
    ax.legend(loc='upper center',bbox_to_anchor=(.5,1.18),frameon=False,fontsize=13,ncol=2)
    save(fig,'plate-comparison.png')

    fluid = source('results/cad_viscous/20261017T201823726/flow.mat')
    with h5py.File(fluid) as f:
        r=f['result']; refs=np.array(r['meanVelocity']).ravel(); u=np.array(f[refs[0]]).T; w=np.array(f[refs[1]]).T
        body=np.array(r['solidMask']).T; uc=(u+np.roll(u,-1,axis=0))/2; wc=(w+np.roll(w,-1,axis=1))/2
        bounds=np.array(r['grid/bounds']).T; spacing=np.array(r['grid/h']).ravel()
    x=np.linspace(bounds[0,0]+spacing[0]/2,bounds[0,1]-spacing[0]/2,u.shape[0])
    z=np.linspace(bounds[1,0]+spacing[1]/2,bounds[1,1]-spacing[1]/2,u.shape[1])
    vort=(np.roll(wc,-1,axis=0)-np.roll(wc,1,axis=0))/(2*spacing[0])-(np.roll(uc,-1,axis=1)-np.roll(uc,1,axis=1))/(2*spacing[1])
    vort[body>.5]=np.nan
    np.savez_compressed(out/'wake-sample.npz',x=x[::2],z=z[::2],vorticity=vort[::2,::2],body=body[::2,::2])
    cmap=LinearSegmentedColormap.from_list('mc-signed',[GREEN,PAPER,RED])
    fig,ax=plt.subplots(figsize=(7.2,4.2));q=ax.pcolormesh(x,z,vort.T,cmap=cmap,vmin=-8,vmax=8,rasterized=True,shading='auto')
    ax.contour(x,z,body.T,levels=[.5],colors=INK,linewidths=.7);ax.set(xlim=(-.25,3),ylim=(-.7,.7),xlabel='x/c',ylabel='z/c');ax.set_aspect('equal')
    fig.colorbar(q,ax=ax,label='Vorticidad transversal ωc/U',extend='both',pad=.025,shrink=.8)
    save(fig,'viscous-wake.png')

    motion = np.loadtxt(out/'posterior-motion.csv',delimiter=','); ref=np.loadtxt(out/'posterior-reference.csv',delimiter=',')
    fig,ax=plt.subplots(figsize=(7.2,4.2));line,=ax.plot(ref[:,0],motion[0,1:],color=GREEN,lw=2)
    ax.set(xlim=(ref[0,0],ref[-1,0]),ylim=(-11,11),xlabel='x/c material · piel superior, η=0,5',ylabel='w respecto al equilibrio [µm]');ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(out/'motion-poster.png',dpi=150)
    plt.rcParams['animation.ffmpeg_path']=imageio_ffmpeg.get_ffmpeg_exe()
    writer=FFMpegWriter(fps=30,codec='libx264',bitrate=900,extra_args=['-pix_fmt','yuv420p','-movflags','+faststart'])
    with writer.saving(fig,str(out/'posterior-motion.mp4'),dpi=120):
        for row in motion:line.set_ydata(row[1:]);writer.grab_frame()
    plt.close(fig)

    snapshot = {'created_date':'2026-10-17','scope':'conditional numerical study; no specimen validation',
                'source_step_sha256':review['sourceStepSha256'],'verification_checks':69,
                'structural':{'dofs':[3*int(r['nodes']) for r in data('results/cad_original/campaign_20261017T190836945/summary.csv')]},
                'speed_map':review['speedMap'],'time_refinement':review['dynamic']['time'],
                'limits':{'curvature_last_change':review['posterior']['curvatureFieldLastRelativeChange'],'viscous_CL_last_change':review['viscous']['lastRelativeCLChange']},
                'sources':sources,'media_license':'MIT for original numerical figures and motion; plate-comparison derived data CC BY 4.0 Huang–Zhang v2',
                'video':{'physical_window_s':[float(motion[0,0]),float(motion[-1,0])],'states':len(motion),'fps':30,'display_duration_s':len(motion)/30,'interpolated_states':False}}
    (out/'evidence.json').write_text(json.dumps(snapshot,indent=2,ensure_ascii=False)+'\n')
    assets=[p for p in out.iterdir() if p.is_file() and p.name!='asset-hashes.json']
    (out/'asset-hashes.json').write_text(json.dumps({p.name:digest(p) for p in sorted(assets)},indent=2)+'\n')
    print('Exported',len(assets),'curated assets to',out)


if __name__=='__main__': main()
