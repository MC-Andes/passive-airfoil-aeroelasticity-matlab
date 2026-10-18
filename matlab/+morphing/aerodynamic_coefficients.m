function r = aerodynamic_coefficients(points,forces,c,U,alpha_deg)
%AERODYNAMIC_COEFFICIENTS Integrate actual point forces about the shaft axis.
% Forces are ON the body in canonical chord/span/height axes, in N.
% Positive drag is along the freestream; positive lift is the wind-normal.
assert(isequal(size(points),size(forces)) && size(points,2)==3 && all(isfinite([points(:);forces(:)])),'morphing:Traction','Require finite physical force points and forces, not Cp alone.');
assert(isscalar(U) && isfinite(U) && U>0 && isscalar(alpha_deg) && isfinite(alpha_deg),'morphing:Flow','Finite positive velocity and incidence required.');
F=sum(forces,1);axis=c.audit.shaft_socket.canonical_center_m;ref=[axis(1),0,axis(3)];M=sum(cross(points-ref,forces,2),1);a=deg2rad(alpha_deg);
q=.5*c.flow.rho_kg_m3*U^2;area=c.geometry.referenceArea_m2;
r=struct('force_N',F,'moment_Nm',M,'reference_m',ref,'q_Pa',q,'Re',c.flow.rho_kg_m3*U*c.geometry.chord_m/c.flow.mu_Pa_s, ...
 'CD',dot(F,[cos(a),0,sin(a)])/(q*area),'CL',dot(F,[-sin(a),0,cos(a)])/(q*area),'Cm',M(2)/(q*area*c.geometry.chord_m), ...
 'convention','freestream [cos(alpha),0,sin(alpha)]; lift [-sin(alpha),0,cos(alpha)]; Cm about +span shaft axis, positive nose-up; half-wing S=.09, c=.25');
end
