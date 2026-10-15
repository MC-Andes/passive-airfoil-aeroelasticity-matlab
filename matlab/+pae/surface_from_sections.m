function [surface,info] = surface_from_sections(folder,nx,ny,edgePolicy)
%SURFACE_FROM_SECTIONS Build mean surface from inspected CAD outer envelopes.
% Canonical axes: x chord, y span, z height. Requires explicit end policy:
% 'constant' extends nearest sampled sections to physical span ends;
% 'sampled' uses only the interior sampled span, with no extrapolation.
arguments
    folder
    nx (1,1) double {mustBeInteger,mustBePositive}
    ny (1,1) double {mustBeInteger,mustBePositive}
    edgePolicy (1,1) string {mustBeMember(edgePolicy,["constant","sampled"])}
end
info=jsondecode(fileread(fullfile(folder,'cad_info.json')));
assert(strcmp(info.units,'m'),'pae:Units','CAD cache must be in metres.');
records=info.sections; n=numel(records); span=zeros(n,1);
sections=cell(n,1);
for j=1:n
    sections{j}=readmatrix(fullfile(folder,records(j).file));
    span(j)=records(j).span_m;
    values=sections{j};
    assert(size(values,2)==5 && all(isfinite(values(:))) && all(diff(values(:,1))>0), ...
        'pae:Section','Invalid sampled section.');
end
assert(n>=2 && all(diff(span)>0),'pae:Section','Invalid section span stations.');
axis=info.cad_axes_zero_based(2)+1;
if edgePolicy=="constant"
    span=[info.minimum_m(axis);span;info.maximum_m(axis)];
    sections=[sections(1);sections;sections(end)];
end
xi=(1-cos(linspace(0,pi,nx+1)))/2;
eta=(1-cos(linspace(0,pi,ny+1)))/2;
targetSpan=span(1)+range(span)*eta;
X=zeros(nx+1,numel(span)); Z=X;
for j=1:numel(span)
    values=sections{j};
    X(:,j)=interp1(values(:,1),values(:,2),xi,'pchip');
    Z(:,j)=interp1(values(:,1),values(:,5),xi,'pchip');
end
nodes=zeros(nx+1,ny+1,3);
nodes(:,:,1)=interp1(span,X',targetSpan,'linear')'-min(X(1,:));
nodes(:,:,2)=repmat(targetSpan-.5*(span(1)+span(end)),nx+1,1);
nodes(:,:,3)=interp1(span,Z',targetSpan,'linear')'-mean(Z(1,:));
[XI,ETA]=ndgrid(xi,eta); chord=X(end,:)-X(1,:);
assert(all(chord>0),'pae:Section','Nonpositive CAD chord.');
surface=struct('nodes',nodes,'xi',XI,'eta',ETA,'nx',nx,'ny',ny, ...
    'referenceArea',trapz(span,chord),'referenceChord',trapz(span,chord.^2)/trapz(span,chord), ...
    'span',range(span),'source','CAD outer-envelope mean surface', ...
    'edgePolicy',edgePolicy,'cadSourceName',info.source_name, ...
    'canonicalOriginCAD',[min(X(1,:)),.5*(span(1)+span(end)),mean(Z(1,:))]);
pae.panels(nodes);
end
