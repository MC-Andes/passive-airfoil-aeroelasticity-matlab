function B = pressure_map(mesh,xEdges,yEdges)
%PRESSURE_MAP Integrate constant pressures on rectangular physical patches.
% Input order: x varies fastest, like pae.panels. Pressure in Pa acts +w.
% Exact overlaps and 2x2 quadrature preserve total force and first moments.
% Avoids a nearest-node footprint that shrinks during structural refinement.
xEdges=xEdges(:); yEdges=yEdges(:);
assert(numel(xEdges)>1 && numel(yEdges)>1 && all(isfinite([xEdges;yEdges])) && ...
    all(diff(xEdges)>0) && all(diff(yEdges)>0),'plate:Pressure','Edges must increase.');
tol=1e-12*max([1;abs(xEdges);abs(yEdges)]);
assert(xEdges(1)>=mesh.x(1)-tol && xEdges(end)<=mesh.x(end)+tol && ...
    yEdges(1)>=mesh.y(1)-tol && yEdges(end)<=mesh.y(end)+tol, ...
    'plate:Pressure','Pressure patches must lie inside structural plate.');
np=(numel(xEdges)-1)*(numel(yEdges)-1); B=sparse(3*size(mesh.nodes,1),np);
g=[-1,1]/sqrt(3);
for e=1:size(mesh.quads,1)
    nodes=mesh.quads(e,:); X=mesh.nodes(nodes,:); lo=min(X); hi=max(X);
    ix=find(xEdges(1:end-1)<hi(1) & xEdges(2:end)>lo(1));
    iy=find(yEdges(1:end-1)<hi(2) & yEdges(2:end)>lo(2));
    for i=ix(:)'
        for j=iy(:)'
            a=max(lo,[xEdges(i),yEdges(j)]); b=min(hi,[xEdges(i+1),yEdges(j+1)]);
            weights=zeros(4,1);
            for u=g
                for v=g
                    point=(a+b)/2+[u,v].*(b-a)/2;
                    xi=2*(point(1)-lo(1))/(hi(1)-lo(1))-1;
                    eta=2*(point(2)-lo(2))/(hi(2)-lo(2))-1;
                    N=plate.shape(xi,eta); weights=weights+N'*prod(b-a)/4;
                end
            end
            panel=i+(j-1)*(numel(xEdges)-1);
            B(3*nodes-2,panel)=B(3*nodes-2,panel)+weights;
        end
    end
end
end
