function [system,fixed] = build_plate_case(c,nChord)
%BUILD_PLATE_CASE Resolve flexible strip and clamp only the rigid root region.
assert(nChord>=4 && nChord==fix(nChord),'ael:Mesh','Use at least four chord cells.');
switchX=c.chord*(1-c.flexibleFraction);
x=unique([linspace(0,c.chord,nChord+1),switchX]);
y=linspace(0,c.span,ceil(c.span/c.chord*nChord)+1);
t=@(x,y) c.rigidThickness+(c.flexibleThickness-c.rigidThickness)*(x>switchX); %#ok<INUSD>
mesh=plate.mesh(x,y,t); system=plate.assemble(mesh,c.material);
fixed=fem3d.dofs(find(abs(mesh.nodes(:,2))<1e-12 & mesh.nodes(:,1)<=switchX+1e-12));
end
