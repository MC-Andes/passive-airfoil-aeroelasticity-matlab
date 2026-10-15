function [faces,owners] = boundary(tetra,nodes)
%BOUNDARY Extract outward-oriented external triangles of a conforming tet mesh.
ne=size(tetra,1);
allfaces=[tetra(:,[2,3,4]);tetra(:,[1,4,3]);tetra(:,[1,2,4]);tetra(:,[1,3,2])];
owner=repmat((1:ne)',4,1);
opposite=[tetra(:,1);tetra(:,2);tetra(:,3);tetra(:,4)];
[~,~,group]=unique(sort(allfaces,2),'rows'); counts=accumarray(group,1);
assert(all(counts<=2),'fem3d:Mesh','Nonmanifold or duplicated tetrahedral faces.');
take=counts(group)==1; faces=allfaces(take,:); owners=owner(take); opposite=opposite(take);
a=nodes(faces(:,1),:); n=cross(nodes(faces(:,2),:)-a,nodes(faces(:,3),:)-a,2);
flip=sum(n.*(nodes(opposite,:)-a),2)>0;
faces(flip,[2,3])=faces(flip,[3,2]);
end
