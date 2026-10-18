function row = response_row(folder,speed)
%RESPONSE_ROW Metrics of one completed speed, including split campaign runs.
subfolder=fullfile(folder,"U_"+speed);s=load(fullfile(subfolder,'dynamic.mat'),'dynamics','metrics','flow');g=jsondecode(fileread(fullfile(subfolder,'growth_fit.json')));
assert(abs(s.flow.speed-speed)<1e-12 && g.qualified,'cad:Review','Speed or growth-fit qualification mismatch.');
row=array2table([speed,g.frequencyHz,g.growthPerSecond,g.peakCount,max(s.metrics.fluctuationRMS),s.dynamics.relativeEnergyBalance,max(s.dynamics.relativeResidual)], ...
 'VariableNames',{'U_ms','frequency_Hz','growth_per_s','peaks','maximum_RMS_m','relative_energy_defect','relative_equilibrium_residual'});
end
