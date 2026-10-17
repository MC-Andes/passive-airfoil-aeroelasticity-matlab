function g = lattice(nodes)
%LATTICE Quarter-chord rings, collocation at 3/4 chord, actual trailing edge.
% Ring order: upstream-left -> upstream-right -> downstream-right -> downstream-left.
% Positive circulation induces downwash and positive lift in positive-x flow.
p=pae.panels(nodes); [nx1,ny1,~]=size(nodes); nx=nx1-1; ny=ny1-1; nn=nx1*ny1;
ids=reshape(1:nn,nx1,ny1); I=[]; J=[]; V=[];
for j=1:ny1
    for i=1:nx
        row=ids(i,j); I=[I;row;row]; J=[J;ids(i,j);ids(i+1,j)]; V=[V;.75;.25]; %#ok<AGROW>
    end
    I=[I;ids(end,j)]; J=[J;ids(end,j)]; V=[V;1]; %#ok<AGROW>
end
G=sparse(I,J,V,nn,nn); R=G*reshape(nodes,[],3);
ring=[reshape(ids(1:end-1,1:end-1),[],1),reshape(ids(1:end-1,2:end),[],1), ...
    reshape(ids(2:end,2:end),[],1),reshape(ids(2:end,1:end-1),[],1)];
spanA=reshape(ids(1:end-1,1:end-1),[],1); spanB=reshape(ids(1:end-1,2:end),[],1);
chordA=reshape(ids(1:end-1,:),[],1); chordB=reshape(ids(2:end,:),[],1);
Hspan=(G(spanA,:)+G(spanB,:))/2; Hchord=(G(chordA,:)+G(chordB,:))/2;
np=nx*ny; rows=repmat((1:np)',1,4); Hcent=sparse(rows(:),p.corners(:),.25,np,nn);
g=struct('nodes',nodes,'ringVertices',R,'ring',ring,'panels',p,'nx',nx,'ny',ny, ...
    'Hspan',Hspan,'Hchord',Hchord,'Hcent',Hcent,'G',G, ...
    'spanA',R(spanA,:),'spanB',R(spanB,:), ...
    'chordA',R(chordA,:),'chordB',R(chordB,:), ...
    'trailingEdge',reshape(nodes(end,:,:),ny1,3));
end
