function value = cad_mask(coordinates,time,profile,span,settings)
%CAD_MASK Closed exterior envelope from an audited CAD section, in chord units.
% 2D coordinates x/z; 3D x/z/span. Internal sealed cavities have no exterior
% fluid boundary here. Open span ends are explicitly capped in the 3D scenario.
x=coordinates{1};z=coordinates{2};alpha=settings.alpha;
angle=-alpha;shift=0; % positive flow incidence corresponds to negative body rotationif isfield(settings,'amplitude'),shift=settings.amplitude*sin(settings.omega*time);end
% Rigid prescribed heave is a wall-motion verification, not flow-excited FSI.
xx=.25+cos(angle)*(x-.25)+sin(angle)*(z-shift);zz=-sin(angle)*(x-.25)+cos(angle)*(z-shift);
upper=interp1(profile(:,1),profile(:,2),xx,'linear',NaN);lower=interp1(profile(:,1),profile(:,3),xx,'linear',NaN);
distance=min(upper-zz,zz-lower);distance(~isfinite(distance))=-inf;
if numel(coordinates)==3,distance=min(distance,span/2-abs(coordinates{3}));end
width=settings.smoothingWidth;value=max(0,min(1,.5+distance/width));
end
