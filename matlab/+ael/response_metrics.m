function metrics = response_metrics(time,Y)
%RESPONSE_METRICS Common trailing-edge observables for a selected time window.
% No automatic transient removal. Dominant FFT bin is not a natural frequency.
% Uniform sampling required; explicit Hann window implemented in MATLAB base.
time=time(:)'; nt=numel(time); dt=diff(time);
assert(nt>=8 && size(Y,2)==nt && all(isfinite([time(:);Y(:)])) && ...
    all(dt>0) && max(abs(dt-mean(dt)))<1e-8*mean(dt), ...
    'ael:Sampling','Require at least eight uniformly sampled finite values.');
average=mean(Y,2); fluctuation=Y-average; w=.5-.5*cos(2*pi*(0:nt-1)/(nt-1));
spectrum=abs(fft(fluctuation.*w,[],2))/sum(w);
frequency=(0:floor(nt/2))/(nt*mean(dt)); spectrum=spectrum(:,1:numel(frequency));
if rem(nt,2)==0, spectrum(:,2:end-1)=2*spectrum(:,2:end-1);
else, spectrum(:,2:end)=2*spectrum(:,2:end); end
[peak,idx]=max(spectrum(:,2:end),[],2); dominant=frequency(idx+1)';
dominant=dominant(:); dominant(max(abs(fluctuation),[],2)<1e-12*max(1,max(abs(Y),[],2)))=NaN;
metrics=struct('mean',average,'fluctuationRMS',sqrt(mean(fluctuation.^2,2)), ...
    'peakToPeak',max(Y,[],2)-min(Y,[],2),'dominantFrequencyHz',dominant, ...
    'frequencyResolutionHz',1/(nt*mean(dt)),'frequencyHz',frequency, ...
    'peakAmplitudeSpectrum',spectrum,'dominantBinPeakAmplitude',peak);
end
