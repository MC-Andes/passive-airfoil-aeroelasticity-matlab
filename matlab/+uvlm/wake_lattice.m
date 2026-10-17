function [vertices,rings] = wake_lattice(nodes)
%WAKE_LATTICE Same ring orientation as the bound surface.
[nr,ny1,~]=size(nodes); vertices=reshape(nodes,[],3);
if nr<2, rings=zeros(0,4); return; end
ids=reshape(1:nr*ny1,nr,ny1);
rings=[reshape(ids(1:end-1,1:end-1),[],1),reshape(ids(1:end-1,2:end),[],1), ...
    reshape(ids(2:end,2:end),[],1),reshape(ids(2:end,1:end-1),[],1)];
end
