function surface = surface_traction(result,profile,span,alpha,offsetCells)
%SURFACE_TRACTION Exterior pressure/viscous stress sampling on an immersed wall.
% Samples offset into fluid to avoid reporting pressure inside the penalty
% layer as resolved skin pressure. Offset is a declared discretization error.
if nargin<5,offsetCells=2;end
g=result.grid;d=g.dimension;assert(d==2,'nsflow:Surface','This section operator is for 2D per-unit-span tractions.');
x=profile(:,1);upper=profile(:,2);lower=profile(:,3);nu=1/result.summary.Re;h=min(g.h);R=[cos(alpha),sin(alpha);-sin(alpha),cos(alpha)];
center=cell(1,2);for a=1:2,center{a}=(result.meanVelocity{a}+nsflow.shift(result.meanVelocity{a},a,-1))/2;end
referenceRegion=g.pressureCoordinates{1}<g.bounds(1,1)+.1*diff(g.bounds(1,:)) & abs(g.pressureCoordinates{2}-mean(g.bounds(2,:)))>.3*diff(g.bounds(2,:));assert(any(referenceRegion,'all'),'nsflow:Reference','No upstream reservoir pressure reference.');referencePressure=mean(result.meanPressure(referenceRegion));
axes={g.pressureCoordinates{1}(:,1),g.pressureCoordinates{2}(1,:)'};P=griddedInterpolant(axes,result.meanPressure,'linear','none');gradient=cell(2,2);
for a=1:2,for b=1:2,derivative=(nsflow.shift(center{a},b,-1)-nsflow.shift(center{a},b,1))/(2*g.h(b));gradient{a,b}=griddedInterpolant(axes,derivative,'linear','none');end,end
rows=[];total=zeros(1,2);pressureForce=total;viscousForce=total;
for side=1:2
 if side==1,z=upper;sign=1;else,z=lower;sign=-1;end
 for k=1:numel(x)-1
  delta=[x(k+1)-x(k),z(k+1)-z(k)];length=norm(delta);if length<eps,continue;end
  normal=sign*[-delta(2),delta(1)]/length;point=[mean(x(k:k+1))-.25,mean(z(k:k+1))]*R'+[.25,0];normal=normal*R';sample=point+offsetCells*h*normal;p=P(sample(1),sample(2))-referencePressure;H=zeros(2);
  for a=1:2,for b=1:2,H(a,b)=gradient{a,b}(sample(1),sample(2));end,end
  pressure=-p*normal;viscous=nu*((H+H')*normal')';force=(pressure+viscous)*length;total=total+force;pressureForce=pressureForce+pressure*length;viscousForce=viscousForce+viscous*length;
  rows(end+1,:)=[side,mean(x(k:k+1)),point,normal,2*p,force,length]; %#ok<AGROW>
 end
end
surface=struct('rows',rows,'totalForcePerSpan',total,'pressureForcePerSpan',pressureForce,'viscousForcePerSpan',viscousForce,'span',span, ...
 'forceCoefficient',2*total,'pressureReferenceValue',referencePressure,'pressureReference','upstream 10% of box, outer 40% of height; reservoir finite-domain reference','offsetCells',offsetCells,'offsetChord',offsetCells*h, ...
 'scope','pressure/shear sampled outside immersed penalty layer; offset and grid convergence required for quantitative skin Cp');
end
