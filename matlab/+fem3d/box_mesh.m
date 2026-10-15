function mesh = box_mesh(lengths,divisions)
%BOX_MESH Rectangular solid, origin at (0,0,0); SI units.
arguments
    lengths (1,3) double {mustBePositive,mustBeFinite}
    divisions (1,3) double {mustBeInteger,mustBePositive}
end
[x,y,z]=ndgrid(linspace(0,lengths(1),divisions(1)+1), ...
    linspace(0,lengths(2),divisions(2)+1),linspace(0,lengths(3),divisions(3)+1));
mesh=fem3d.structured_mesh(cat(4,x,y,z));
end
