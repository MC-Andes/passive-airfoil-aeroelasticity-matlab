function result = harmonic(model,frequencyHz,input)
%HARMONIC Complex amplitudes for exp(i*omega*t); prescribed pressure/force only.
f=frequencyHz(:)';
assert(all(isfinite(f) & f>=0) && size(input,1)==size(model.B,2) && ...
    (size(input,2)==1 || size(input,2)==numel(f)) && all(isfinite(input(:))), ...
    'ael:Frequency','Invalid frequency or input amplitudes.');
if size(input,2)==1, input=repmat(input,1,numel(f)); end
Q=complex(zeros(size(model.K,1),numel(f)));
for k=1:numel(f)
    w=2*pi*f(k); Z=model.K-w^2*model.M+1i*w*model.C;
    assert(rcond(full(Z))>1e-14,'ael:Resonance','Singular undamped resonance.');
    Q(:,k)=Z\(model.B*input(:,k));
end
result=struct('frequencyHz',f,'Q',Q,'Y',model.O*Q, ...
    'convention','Re(Y*exp(i*2*pi*f*t)); amplitudes are peak, not RMS');
end
