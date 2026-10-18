function V = segment_influence(points,A,B,core)
%SEGMENT_INFLUENCE Unit Biot-Savart for each endpoint pair, without summation.
% Output Npoints x 3 x Nsegments; bounded batches are chosen by callers.
if isreal(points) && isreal(A) && isreal(B) && exist('aero_vortex_kernel','file')==3
    V=aero_vortex_kernel(0,points,A,B,core);return;
end
np=size(points,1);ns=size(A,1);r0=reshape(B-A,1,ns,3);r1=reshape(points,np,1,3)-reshape(A,1,ns,3);r2=reshape(points,np,1,3)-reshape(B,1,ns,3);
n1=sqrt(sum(r1.^2,3));n2=sqrt(sum(r2.^2,3));product=cross(r1,r2,3);den=sum(product.^2,3)+core^2*sum(r0.^2,3);valid=real(n1)>core*1e-3 & real(n2)>core*1e-3 & real(den)>realmin;
n1(~valid)=1;n2(~valid)=1;den(~valid)=1;factor=sum((r1./n1-r2./n2).*r0,3)./den;factor(~valid)=0;V=permute(product.*factor/(4*pi),[1,3,2]);
end
