function D = isotropic(E,nu)
%ISOTROPIC Linear elastic matrix, SI; engineering shear strains.
% Order: [xx yy zz xy yz xz]; stress shear components are not doubled.
arguments
    E (1,1) double {mustBePositive,mustBeFinite}
    nu (1,1) double {mustBeFinite}
end
assert(nu>-1 && nu<.5,'fem3d:Material','Require -1 < nu < 0.5.');
lambda=E*nu/((1+nu)*(1-2*nu)); mu=E/(2*(1+nu));
D=zeros(6); D(1:3,1:3)=lambda; D(1:3,1:3)=D(1:3,1:3)+2*mu*eye(3);
D(4:6,4:6)=mu*eye(3);
end
