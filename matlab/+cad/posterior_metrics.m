function metrics = posterior_metrics(points,displacement,xi,eta)
%POSTERIOR_METRICS Material upper-skin shape, x-height slope/curvature change.
% Cubic fits over a fixed +/-0.025 chord window avoid differentiating mesh
% facets or the singular geometric trailing edge. Curvature units are 1/m.
ny=numel(eta);nx=numel(xi);assert(size(points,1)==ny*nx && isequal(size(points),size(displacement)),'cad:Shape','Material sampling dimensions differ.');
reference=reshape(points,ny,nx,3);deformed=reshape(points+displacement,ny,nx,3);centers=find(xi>=xi(1)+.025-1e-12 & xi<=xi(end)-.025+1e-12);rows=[];slope=zeros(ny,numel(centers));curvature=slope;
for i=1:numel(centers)
 center=centers(i);use=abs(xi-xi(center))<=.025+1e-12;d=xi(use)'-xi(center);A=[ones(numel(d),1),d,d.^2,d.^3];
 for j=1:ny
  a=A\reshape(reference(j,use,:),[],3);b=A\reshape(deformed(j,use,:),[],3);ka=(a(2,1)*2*a(3,3)-a(2,3)*2*a(3,1))/(a(2,1)^2+a(2,3)^2)^(3/2);kb=(b(2,1)*2*b(3,3)-b(2,3)*2*b(3,1))/(b(2,1)^2+b(2,3)^2)^(3/2);
  slope(j,i)=atan2(b(2,3),b(2,1))-atan2(a(2,3),a(2,1));curvature(j,i)=kb-ka;rows(end+1,:)=[xi(center),eta(j),slope(j,i),curvature(j,i),ka,kb]; %#ok<AGROW>
 end
end
metrics=struct('slopeChange_rad',slope,'curvatureChange_per_m',curvature,'centers',xi(centers),'rows',rows,'fitHalfWindow_over_c',.025,'definition','signed curvature of material-line projection onto canonical chord/height plane; cubic physical window, not element-facet derivatives');
end
