function V = ring_velocity(points,vertices,rings,gamma,core)
%RING_VELOCITY Vectorized Biot-Savart sum with a finite vortex core.
% Sources are processed in blocks. No circulation is deleted or decayed.
P=points; np=size(P,1); V=zeros(np,3); gamma=gamma(:);
assert(numel(gamma)==size(rings,1) && all(isfinite(gamma)),'uvlm:Circulation','Invalid ring circulation.');
keep=gamma~=0; rings=rings(keep,:); gamma=gamma(keep);
if isempty(gamma), return; end
for start=1:128:size(rings,1)
    ix=start:min(start+127,size(rings,1)); ri=rings(ix,:); gam=reshape(gamma(ix),1,[],1);
    for edge=1:4
        A=vertices(ri(:,edge),:); B=vertices(ri(:,mod(edge,4)+1),:);
        if exist('aero_vortex_kernel','file')==3
            V=V+sum(uvlm.segment_influence(P,A,B,core).*reshape(gamma(ix),1,1,[]),3);continue;
        end
        r0=reshape(B-A,1,[],3); r1=reshape(P,np,1,3)-reshape(A,1,[],3); r2=reshape(P,np,1,3)-reshape(B,1,[],3);
        n1=sqrt(sum(r1.^2,3)); n2=sqrt(sum(r2.^2,3)); cr=cross(r1,r2,3);
        den=sum(cr.^2,3)+core^2*sum(r0.^2,3);
        valid=n1>core*1e-3 & n2>core*1e-3 & den>realmin;
        factor=sum((r1./max(n1,realmin)-r2./max(n2,realmin)).*r0,3)./max(den,realmin);
        factor(~valid)=0;
        V=V+reshape(sum(cr.*factor.*gam,2),np,3)/(4*pi);
    end
end
end
