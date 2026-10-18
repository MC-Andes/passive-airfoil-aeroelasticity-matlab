function report = resource_estimate(mesh)
%RESOURCE_ESTIMATE Count actual edges before allocating a quadratic CAD system.
T=mesh.tetra(:,1:4); ne=size(T,1); edges=[1,2;2,3;3,1;1,4;2,4;3,4]; all=zeros(6*ne,2);
for k=1:6, all((k-1)*ne+(1:ne),:)=sort(T(:,edges(k,:)),2); end
nedge=size(unique(all,'rows'),1); nn=size(mesh.nodes,1); nodes10=nn+nedge;
report=struct('tetrahedra',ne,'vertices',nn,'uniqueEdges',nedge,'tet10Nodes',nodes10,'tet10Dofs',3*nodes10, ...
    'uncompressedKMTripletUpperBytes',ne*900*32,'assemblyChunkTripletBytes',floor(432000/900)*900*32, ...
    'factorizationBudget','not estimated by element count; pilot factor nonzeros and process RSS must be measured', ...
    'geometry','straight edge elevation does not improve CAD curvature or sliver quality');
end
