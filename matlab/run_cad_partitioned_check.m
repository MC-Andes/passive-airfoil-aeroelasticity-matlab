function report = run_cad_partitioned_check(bundlePath)
%RUN_CAD_PARTITIONED_CHECK Actual CAD full UVLM strong-FSI vs tangent coupling.
root=fileparts(fileparts(mfilename('fullpath')));s=load(bundlePath,'bundle');b=s.bundle;m=b.model;flow=b.flow;uv=m;uv.aero=@fsi.uvlm_static;eq=fsi.static(uv,flow);m.X0=reshape(reshape(m.X0,[],3)+reshape(m.PhiSurface*eq.q,3,[])',size(m.X0));
time=0:.0005:.1;options=struct('wakeLength',1000*flow.referenceChord);fluid=uvlm.linear_model(m.X0,m.PhiSurface,flow,.0005,numel(time),options);q0=zeros(size(m.K,1),1);q0(1)=1e-5/max(abs(m.observation(:,1)));v0=zeros(size(q0));
g=fluid.state0.geometry;N=zeros(numel(m.X0),size(m.K,1));for j=1:size(m.K,1),nodal=g.Hcent'*(flow.rho*g.panels.area.*g.panels.normal.*(fluid.E(:,j)));N(:,j)=reshape(nodal',[],1);end
m.nodalAccelerationLoad=N;m.initializeFluid=@initialize;m.advanceFluid=@advance;baseline=fluid.base.nodalForces;
linear=fsi.linear_dynamic(m,fluid,time,q0,v0);a=fsi.dynamic(m,flow,time,q0,v0,struct('relaxation',.6,'relativeTolerance',1e-6,'forceScale',norm(m.K*q0)));bb=fsi.dynamic(m,flow,time,q0,v0,struct('relaxation',.8,'relativeTolerance',1e-8,'forceScale',norm(m.K*q0)));Y=m.observation*linear.Q;Ya=m.observation*a.Q;Yb=m.observation*bb.Q;
report=struct('nonlinearTangentSensorRelative',max(abs(Yb-Y),[],'all')/1e-5,'relaxationToleranceSensorRelative',max(abs(Ya-Yb),[],'all')/1e-5,'energyDefect',bb.relativeEnergyBalance,'interfaceWorkDefect_J',bb.maximumInterfaceWorkDefect,'maximumSubiterations',max(bb.iterations),'steps',numel(time)-1,'scope','actual CAD 8x6 material surface, full nonlinear UVLM with strong subiterations vs its tangent; prescribed initial perturbation, no harmonic forcing');
output=fullfile(root,'results/cad_original/partitioned_check');if ~isfolder(output),mkdir(output);end;cad.write_json(fullfile(output,'comparison.json'),report);save(fullfile(output,'responses.mat'),'linear','a','bb','Y','Ya','Yb','report','-v7.3');disp(report);assert(report.nonlinearTangentSensorRelative<.02 && report.relaxationToleranceSensorRelative<.02 && report.energyDefect<.02 && report.interfaceWorkDefect_J<1e-10);
 function [state,aero]=initialize(X,V,f),[state,aero]=uvlm.initialize(X,V,f,options);aero.nodalForces=aero.nodalForces-baseline;end
 function [state,aero]=advance(old,X,V,f,dt),[state,aero]=uvlm.step(old,X,V,f,dt);aero.nodalForces=aero.nodalForces-baseline;end
end
