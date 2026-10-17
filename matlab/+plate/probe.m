function O = probe(mesh,points,component)
%PROBE Bilinear observation at fixed physical points, without nearest-node drift.
% points [x y] metres; component 1=w, 2=theta_x, 3=theta_y.
assert(size(points,2)==2 && all(isfinite(points(:))) && ...
    isscalar(component) && ismember(component,1:3),'plate:Probe','Invalid probes.');
O=sparse(size(points,1),3*size(mesh.nodes,1));
for k=1:size(points,1)
    p=points(k,:); x=mesh.x; y=mesh.y;
    tol=1e-12*max([1;abs(x);abs(y)]);
    assert(p(1)>=x(1)-tol && p(1)<=x(end)+tol && p(2)>=y(1)-tol && ...
        p(2)<=y(end)+tol,'plate:Probe','Probe outside plate.');
    i=find(x<=p(1)+tol,1,'last'); i=min(i,numel(x)-1);
    j=find(y<=p(2)+tol,1,'last'); j=min(j,numel(y)-1);
    xi=2*(p(1)-x(i))/(x(i+1)-x(i))-1; eta=2*(p(2)-y(j))/(y(j+1)-y(j))-1;
    N=plate.shape(xi,eta); e=i+(j-1)*(numel(x)-1);
    O(k,3*mesh.quads(e,:)-3+component)=N;
end
end
