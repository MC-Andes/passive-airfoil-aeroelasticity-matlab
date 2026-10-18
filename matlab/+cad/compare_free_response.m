function comparison=compare_free_response(coarseFolder,fineFolder,speed)
if nargin<3,speed=10;end;subfolder="U_"+speed;
a=load(fullfile(coarseFolder,subfolder,'dynamic.mat'),'dynamics','Y');b=load(fullfile(fineFolder,subfolder,'dynamic.mat'),'dynamics','Y');
ga=jsondecode(fileread(fullfile(coarseFolder,subfolder,'growth_fit.json')));gb=jsondecode(fileread(fullfile(fineFolder,subfolder,'growth_fit.json')));
t=a.dynamics.time;Y=interp1(b.dynamics.time',b.Y',t','linear')';q=interp1(b.dynamics.time',b.dynamics.Q(1,:)',t','linear')';
history=norm(a.Y-Y,inf)/max(norm(Y,inf),realmin); % induced infinity norm, time-integrated row sum
peak=max(abs(a.Y-Y),[],'all')/max(abs(Y),[],'all');
ra=sqrt(mean(a.Y.^2,2));rb=sqrt(mean(Y.^2,2));rmsChange=norm(ra-rb,inf)/max(norm(rb,inf),realmin);
use=t>=.15 & t<=.65;basis=exp(gb.growthPerSecond*t(use)').*[cos(2*pi*gb.frequencyHz*t(use)'),sin(2*pi*gb.frequencyHz*t(use)')];
ca=basis\a.dynamics.Q(1,use)';cb=basis\q(use)';za=ca(1)-1i*ca(2);zb=cb(1)-1i*cb(2);phase=abs(rad2deg(angle(za/zb)));amplitude=abs(abs(za/zb)-1);
frequency=abs(ga.frequencyHz/gb.frequencyHz-1);growth=abs(ga.growthPerSecond/gb.growthPerSecond-1);
comparison=struct('coarse',coarseFolder,'fine',fineFolder,'frequencyRelative',frequency,'growthRelative',growth,'maximumSensorHistoryRelative',peak, ...
 'sensorRowSumRelative',history,'rmsRelative',rmsChange,'fittedAmplitudeRelative',amplitude,'phase_deg',phase,'phaseWindow_s',[min(t(use)),max(t(use))], ...
 'passed',frequency<.01 && growth<.02 && peak<.02 && rmsChange<.02 && phase<2);
end
