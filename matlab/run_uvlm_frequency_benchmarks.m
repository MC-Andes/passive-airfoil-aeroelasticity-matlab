function report = run_uvlm_frequency_benchmarks(k)
%RUN_UVLM_FREQUENCY_BENCHMARKS Exact periodic solution of discrete own UVLM.
% Independent Theodorsen heave and Sears convected-gust formulas. The wake
% influence is built by the same physical Biot-Savart kernel; no fitted loads.
if nargin<1,k=.5;end
root=fileparts(fileparts(mfilename('fullpath')));rows=[];omega=2*k;amplitude=.001;aspect=1000;
C=besselh(1,2,k)/(besselh(1,2,k)+1i*besselh(0,2,k));S=(besselj(0,k)-1i*besselj(1,k))*C+1i*besselj(1,k);
expected=[2*pi*k^2*amplitude-4*pi*1i*k*C*amplitude,2*pi*amplitude*S*exp(-1i*k)];
settings=[16,1/64,100,1e-5;32,1/128,100,1e-5;64,1/256,100,1e-5;128,1/512,100,1e-5;64,1/512,100,1e-5;64,1/512,200,1e-5;64,1/512,100,1e-6;128,1/1024,100,1e-5;128,1/1024,200,1e-5;128,1/1024,100,1e-6];
if k~= .5,settings=[128,1/512,100,1e-5;256,1/1024,100,1e-5;256,1/2048,100,1e-5;256,1/2048,200,1e-5;256,1/2048,100,1e-6];end
if k>3,settings=[256,1/1024,100,1e-5;512,1/2048,100,1e-5;512,1/4096,100,1e-5;512,1/4096,200,1e-5;512,1/4096,100,1e-6];end
for test=1:size(settings,1)
 nx=settings(test,1);dt=settings(test,2);horizon=settings(test,3);core=settings(test,4);
 p=struct('spanStations',[-aspect/2;aspect/2],'leadingEdge',[0,-aspect/2,0;0,aspect/2,0],'chord',[1;1]);s=pae.make_surface(p,nx,1);s.nodes(:,:,1)=repmat(linspace(0,1,nx+1)',1,2);g=uvlm.lattice(s.nodes);A=uvlm.influence(g,core);
 age=(0:ceil(horizon/dt))';N=numel(age)-1;wake=reshape(g.trailingEdge,1,2,3)+reshape(age*dt,[],1,1).*reshape([1,0,0],1,1,3);
 [vertices,rings]=uvlm.wake_lattice(wake);gamma=exp(-1i*omega*dt*(1:N)');V=uvlm.ring_velocity(g.panels.control,vertices,rings,gamma,core);W=sum(V.*g.panels.normal,2);TE=sparse(1,nx,1,1,nx);
 rhs=[repmat(1i*omega*amplitude,nx,1),-amplitude*exp(-1i*omega*g.panels.control(:,1))];G=(A+W*TE)\rhs;
 span=G-[zeros(1,2);G(1:end-1,:)];derivative=(1-exp(-1i*omega*dt))/dt;
 CL=2*(sum(span,1)*aspect+sum(g.panels.area.*G,1)*derivative)/aspect;
 for problem=1:2
  error=abs(CL(problem)/expected(problem)-1);ampError=abs(abs(CL(problem)/expected(problem))-1);phaseError=abs(rad2deg(angle(CL(problem)/expected(problem))));
  rows(end+1,:)=[problem,nx,dt,horizon,core,real(CL(problem)),imag(CL(problem)),error,ampError,phaseError]; %#ok<AGROW>
 end
 fprintf('UVLM periodic n=%d dt=%g: heave %.3g, gust %.3g complex error\n',nx,dt,rows(end-1,8),rows(end,8));
end
report=array2table(rows,'VariableNames',{'problem_heave1_gust2','nx','dt','wake_chords','core_over_c','CL_real','CL_imag','complex_relative_error','amplitude_relative_error','phase_error_deg'});output=fullfile(root,'results/cad_original');writetable(report,fullfile(output,sprintf('uvlm-frequency-k%g.csv',k)));
assert(all(report.amplitude_relative_error(end-3:end)<.02) && all(report.phase_error_deg(end-3:end)<2),'uvlm:Benchmark','Refined frequency response exceeds the declared amplitude/phase budget.');
end
