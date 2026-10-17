function [faces,owners] = boundary(tetra,nodes)
%BOUNDARY Outward Tri3/Tri6 faces; Tet10 edges [12 23 31 14 24 34].
ne=size(tetra,1); order=size(tetra,2);
assert(ismember(order,[4,10]),'fem3d:Mesh','Require Tet4 or Tet10 connectivity.');
local=[2,3,4;1,4,3;1,2,4;1,3,2];
if order==10, local=[local,[6,10,9;8,10,7;5,9,8;7,6,5]]; end
allfaces=zeros(4*ne,size(local,2));
for k=1:4, allfaces((k-1)*ne+(1:ne),:)=tetra(:,local(k,:)); end
owner=repmat((1:ne)',4,1); opposite=reshape(tetra(:,1:4),[],1);
[~,~,group]=unique(sort(allfaces(:,1:3),2),'rows'); counts=accumarray(group,1);
assert(all(counts<=2),'fem3d:Mesh','Nonmanifold or duplicated tetrahedral faces.');
if order==10
    [~,~,fullGroup]=unique(sort(allfaces,2),'rows'); fullCounts=accumarray(fullGroup,1);
    assert(all(fullCounts(fullGroup)==counts(group)),'fem3d:Mesh','Nonconforming quadratic face nodes.');
end
take=counts(group)==1; faces=allfaces(take,:); owners=owner(take); opposite=opposite(take);
a=nodes(faces(:,1),:); n=cross(nodes(faces(:,2),:)-a,nodes(faces(:,3),:)-a,2);
flip=sum(n.*(nodes(opposite,:)-a),2)>0;
if order==4, perm=[1,3,2]; else, perm=[1,3,2,6,5,4]; end
faces(flip,:)=faces(flip,perm);
end
