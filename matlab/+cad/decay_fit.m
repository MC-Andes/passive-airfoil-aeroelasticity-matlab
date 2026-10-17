function fit = decay_fit(time,Y)
%DECAY_FIT Positive peaks of a free response; no synthetic flow excitation.
% A finite record gives a growth estimate, never an automatic flutter claim.
time=time(:)';y=Y(:)';indices=find(y(2:end-1)>y(1:end-2) & y(2:end-1)>=y(3:end))+1;
if numel(indices)<4,fit=struct('frequencyHz',NaN,'growthPerSecond',NaN,'peakCount',numel(indices),'qualified',false);return;end
amplitude=abs(y(indices));good=amplitude>max(amplitude)*1e-3;indices=indices(good);amplitude=amplitude(good);
% Exclude two initial peaks, then stop before the fundamental decays into
% faster-mode contamination or numerical noise. Keep the window explicit.
if numel(indices)>5,indices=indices(3:end);amplitude=amplitude(3:end);end
t=time(indices);A=[ones(numel(t),1),t'];coefficient=A\log(amplitude(:));r=log(amplitude(:))-A*coefficient;
fit=struct('frequencyHz',1/mean(diff(t)),'growthPerSecond',coefficient(2),'peakCount',numel(indices), ...
    'amplitudeCutoffRelative',1e-3,'fitWindow_s',[t(1),t(end)],'peakTimes',t,'peakAmplitude',amplitude,'logEnvelopeResidualRMS',sqrt(mean(r.^2)), ...
    'qualified',numel(indices)>=15,'interpretation','growth/decay from free perturbation; check time, panels, wake, modes and damping before instability attribution');
end
