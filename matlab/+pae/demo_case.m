function [model,flow] = demo_case(nx,ny)
%DEMO_CASE Parametric benchmark with demonstrative passive compliance.
if nargin<1, nx=4; end
if nargin<2, ny=12; end
planform=struct('spanStations',[-.18;.18],'leadingEdge',[0,-.18,0;0,.18,0], ...
    'chord',[.25;.25],'camberFcn',@(xi,eta) .02*4*xi.*(1-xi)+0*eta);
surface=pae.make_surface(planform,nx,ny);
model=pae.make_model(surface,pae.example_structure(surface));
flow=struct('rho',1.225,'speed',12,'alpha',deg2rad(5),'beta',0, ...
    'referenceArea',surface.referenceArea,'referenceChord',surface.referenceChord, ...
    'referencePoint',[surface.referenceChord/4,0,0]);
end
