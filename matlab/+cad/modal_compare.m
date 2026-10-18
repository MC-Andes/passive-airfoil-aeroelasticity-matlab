function report = modal_compare(a,b,eta)
%MODAL_COMPARE MAC in common material sensors, all three translations.
% Weight each station by its represented span interval. No FE-node matching.
% Close/repeated modes are flagged; their individual pairing is not a closure test.
w=diff([eta(1);(eta(1:end-1)+eta(2:end))/2;eta(end)]); W=repelem(w,3);
A=a.observed; B=b.observed; crossMAC=abs(A'*(W.*B)).^2./ ...
    (sum(W.*abs(A).^2,1)'*sum(W.*abs(B).^2,1));
n=min(size(A,2),size(B,2)); assert(n<=8,'cad:Modes','Assignment enumerates at most eight modes.');
p=perms(1:n); score=zeros(size(p,1),1);
for k=1:size(p,1), score(k)=sum(crossMAC(sub2ind([n,n],1:n,p(k,:)))); end
[~,best]=max(score); pairing=p(best,:)'; mac=crossMAC(sub2ind([n,n],(1:n)',pairing));
relative=abs(b.frequencyHz(pairing)./a.frequencyHz(1:n)-1);
near=false(n,1);
for k=1:n
    other=setdiff(1:n,k);
    near(k)=any(abs(a.frequencyHz(other)/a.frequencyHz(k)-1)<.01) || ...
        any(abs(b.frequencyHz(setdiff(1:n,pairing(k)))/b.frequencyHz(pairing(k))-1)<.01);
end
report=struct('pairing',pairing,'MAC',mac,'frequencyRelativeChange',relative,'closeModes',near, ...
    'MACmatrix',crossMAC,'weightDefinition','trapezoid span weights on xyz at nine common material stations', ...
    'requiresSubspaceComparison',any(near),'individualPairingReliable',all(mac>.9) && ~any(near));
end
