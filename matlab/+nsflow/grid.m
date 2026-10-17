function g = grid(count,bounds)
%GRID Periodic Cartesian MAC layout in two or three dimensions.
% Pressure at centers; velocity i at the lower face in direction i.
count=count(:)';d=numel(count);assert(ismember(d,[2,3]) && isequal(size(bounds),[d,2]) && ...
    all(count>=8 & count==fix(count)) && all(diff(bounds,1,2)>0),'nsflow:Grid','Invalid periodic Cartesian grid.');
h=(bounds(:,2)-bounds(:,1))'./count;centers=cell(1,d);faces=cell(d,d);
for j=1:d
    centers{j}=bounds(j,1)+((0:count(j)-1)+.5)*h(j);
    for i=1:d,faces{i,j}=centers{j};if i==j,faces{i,j}=faces{i,j}-.5*h(j);end,end
end
P=cell(1,d);[P{:}]=ndgrid(centers{:});X=cell(d,d);
for i=1:d,[X{i,:}]=ndgrid(faces{i,:});end
lambda=zeros(count);
for j=1:d
    shape=ones(1,d);shape(j)=count(j);lambda=lambda+reshape(-4*sin(pi*(0:count(j)-1)/count(j)).^2/h(j)^2,shape);
end
lambda(1)=1;
g=struct('count',count,'dimension',d,'bounds',bounds,'h',h,'pressureCoordinates',{P}, ...
    'velocityCoordinates',{X},'laplacianEigenvalues',lambda,'cellVolume',prod(h));
end
