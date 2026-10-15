function mesh = airfoil_mesh(chord,span,thicknessRatio,divisions)
%AIRFOIL_MESH Solid NACA 00xx surrogate, not the internal reference mechanism.
% Cut at x/c=.005 and .995 to avoid zero-volume leading/trailing cells.
% x chord, y span, z height; centred span. No cavities or slider.
arguments
    chord (1,1) double {mustBePositive,mustBeFinite}
    span (1,1) double {mustBePositive,mustBeFinite}
    thicknessRatio (1,1) double {mustBePositive,mustBeFinite}
    divisions (1,3) double {mustBeInteger,mustBePositive} = [30,16,4]
end
xi=linspace(.005,.995,divisions(1)+1);
[s,y,eta]=ndgrid(xi,linspace(-span/2,span/2,divisions(2)+1), ...
    linspace(-1,1,divisions(3)+1));
half=5*thicknessRatio*chord*(.2969*sqrt(s)-.1260*s-.3516*s.^2+.2843*s.^3-.1015*s.^4);
mesh=fem3d.structured_mesh(cat(4,chord*s,y,eta.*half));
mesh.source="solid truncated NACA 00xx surrogate; no reference internal mechanism";
mesh.chord=chord; mesh.span=span;
end
