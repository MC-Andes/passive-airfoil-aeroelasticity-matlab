function mesh = mesh(x,y,thickness)
%MESH Rectangular flat plate; nonuniform grid must resolve thickness jumps.
% thickness: positive scalar or @(xCentroid,yCentroid), one value per element.
x=x(:); y=y(:);
assert(numel(x)>1 && numel(y)>1 && all(isfinite([x;y])) && ...
    all(diff(x)>0) && all(diff(y)>0),'plate:Mesh','Grid vectors must increase.');
[X,Y]=ndgrid(x,y); ids=reshape(1:numel(X),size(X));
Q=[reshape(ids(1:end-1,1:end-1),[],1),reshape(ids(2:end,1:end-1),[],1), ...
    reshape(ids(2:end,2:end),[],1),reshape(ids(1:end-1,2:end),[],1)];
nodes=[X(:),Y(:)]; centers=zeros(size(Q,1),2);
for k=1:4, centers=centers+nodes(Q(:,k),:)/4; end
if isa(thickness,'function_handle')
    t=thickness(centers(:,1),centers(:,2));
else
    t=thickness;
end
if isscalar(t), t=repmat(t,size(Q,1),1); end
assert(numel(t)==size(Q,1) && all(isfinite(t(:)) & t(:)>0), ...
    'plate:Thickness','Thickness must be positive per element.');
mesh=struct('nodes',nodes,'quads',Q,'thickness',t(:),'centers',centers,'x',x,'y',y);
end
