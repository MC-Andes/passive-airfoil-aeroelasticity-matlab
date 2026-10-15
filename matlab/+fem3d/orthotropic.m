function D = orthotropic(E,nu,G)
%ORTHOTROPIC Constitutive matrix in GLOBAL material axes (no auto rotation).
% E=[Ex Ey Ez], nu=[nu_xy nu_yz nu_xz], G=[Gxy Gyz Gxz].
% nu_ij means lateral strain in j under uniaxial loading in i.
% Reciprocal Poisson ratios follow symmetry of compliance.
arguments
    E (1,3) double {mustBePositive,mustBeFinite}
    nu (1,3) double {mustBeFinite}
    G (1,3) double {mustBePositive,mustBeFinite}
end
S=diag([1./E,1./G]);
S(1,2)=-nu(1)/E(1); S(2,1)=S(1,2);
S(2,3)=-nu(2)/E(2); S(3,2)=S(2,3);
S(1,3)=-nu(3)/E(1); S(3,1)=S(1,3);
[~,p]=chol(S);
assert(p==0,'fem3d:Material','Orthotropic compliance must be positive definite.');
D=S\eye(6); D=(D+D')/2;
end
