function structure = example_structure(surface)
%EXAMPLE_STRUCTURE ILLUSTRATIVE parameters, not identified from STEP.
% Coordinates are metres: camber, internal slider, twist proxy, bending.
% Slider has no direct aerodynamic displacement. K couples it to camber.
N=numel(surface.nodes)/3;
xi=surface.xi; eta=surface.eta;
camber=max(0,(xi-.4)/.6).^2;
spanEnvelope=(2*eta-1).^2;
basis=zeros(N,3,4);
basis(:,3,1)=camber(:);
% Column 2 describes an internal coordinate and therefore has zero skin motion.
basis(:,3,3)=reshape((xi-.25).*spanEnvelope,[],1);
basis(:,3,4)=spanEnvelope(:);
M=diag([.02,.01,.035,.08]);
kCoupling=800; kSkin=1500; kSlider=250;
K=[kSkin+kCoupling,-kCoupling,0,0;-kCoupling,kCoupling+kSlider,0,0; ...
    0,0,6000,0;0,0,0,10000];
% Rayleigh coefficients give 3% damping at the lowest and highest dry modes.
omega=sqrt(sort(eig(K,M))); zeta=.03;
rayleigh=[1/(2*omega(1)),omega(1)/2;1/(2*omega(end)),omega(end)/2]\[zeta;zeta];
C=rayleigh(1)*M+rayleigh(2)*K;
structure=struct('M',M,'C',C,'K',K,'Phi',reshape(basis,3*N,4), ...
    'cubic',[1e8;2e7;1e8;1e8], ...
    'labels',{{'camber','slider','twist_proxy','bending'}}, ...
    'coordinateUnits',{{'m','m','m','m'}}, ...
    'source','ILUSTRATIVO: parámetros sin calibrar; no identificados a partir del CAD');
end
