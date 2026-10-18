function derivative = wake_vertex_gradient(points,centers,segment,core)
%WAKE_VERTEX_GRADIENT Exact endpoint derivative of two collinear vortex legs.
% The shared material wake vertex moves; its two neighbors stay fixed. Same
% regularized Biot-Savart law as segment_influence. Reuses real kernel terms
% for all three coordinate derivatives, avoiding repeated complex arrays.
if exist('aero_vortex_kernel','file')==3
    derivative=cell(1,3);[derivative{:}]=aero_vortex_kernel(1,points,centers,segment,core);return;
end
np=size(points,1);nb=size(centers,1);d=reshape(segment,1,1,3);r=reshape(points,np,1,3)-reshape(centers,1,nb,3);rp=r+d;rm=r-d;
n0=sqrt(sum(r.^2,3));nplus=sqrt(sum(rp.^2,3));nminus=sqrt(sum(rm.^2,3));n0=max(n0,realmin);nplus=max(nplus,realmin);nminus=max(nminus,realmin);
c=cat(3,d(:,:,2).*r(:,:,3)-d(:,:,3).*r(:,:,2),d(:,:,3).*r(:,:,1)-d(:,:,1).*r(:,:,3),d(:,:,1).*r(:,:,2)-d(:,:,2).*r(:,:,1));D=sum(c.^2,3)+core^2*sum(d.^2,3);valid=D>realmin & n0>core*1e-3 & nplus>core*1e-3 & nminus>core*1e-3;D=max(D,realmin);
Splus=sum(d.*(rp./nplus-r./n0),3);Sminus=sum(d.*(r./n0-rm./nminus),3);R=-Splus.*rp+Sminus.*rm;
G=(rp./nplus-2*r./n0+rm./nminus)./D-(2*cross(c,R,3)+2*core^2*d.*(Splus-Sminus))./D.^2;
zero=zeros(np,nb);rotation={cat(3,zero,R(:,:,3),-R(:,:,2)),cat(3,-R(:,:,3),zero,R(:,:,1)),cat(3,R(:,:,2),-R(:,:,1),zero)};derivative=cell(1,3);
for axis=1:3,V=(rotation{axis}./D+c.*G(:,:,axis))/(4*pi);V(repmat(~valid,1,1,3))=0;derivative{axis}=permute(V,[1,3,2]);end
end
