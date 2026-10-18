function [K,M,P,area] = mitc4(X,t,E,nu,rho)
%MITC4 Flat Reissner-Mindlin plate with tied covariant shear strains.
% DOFs per node: [w, theta_x, theta_y]; theta are slopes, gamma=grad(w)-theta.
% X 4x2 in metres; t [m], E [Pa], rho [kg/m^3]. Positive pressure loads +w.
% Bathe & Dvorkin (1985), DOI 10.1002/nme.1620210213. Small rotations only.
assert(isequal(size(X),[4,2]) && all(isfinite(X(:))) && ...
    all(isfinite([t,E,nu,rho])) && t>0 && E>0 && rho>0 && nu>-1 && nu<.5, ...
    'plate:Input','Invalid geometry or isotropic material.');
Db=E*t^3/(12*(1-nu^2))*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
Ds=(5/6)*E/(2*(1+nu))*t*eye(2);
% Tying points: xi-shear on bottom/top, eta-shear on left/right edges.
bottom=tied(0,-1,1,X); top=tied(0,1,1,X);
left=tied(-1,0,2,X); right=tied(1,0,2,X);
K=zeros(12); M=K; P=zeros(12,1); area=0;
g=[-1,1]/sqrt(3);
for xi=g
    for eta=g
        [N,dN]=plate.shape(xi,eta); J=dN*X; determinant=det(J);
        assert(determinant>1e-14*max(vecnorm(X-X(1,:),2,2))^2, ...
            'plate:Mesh','Quad must be nondegenerate and counterclockwise.');
        G=J\dN; B=zeros(3,12);
        B(1,2:3:end)=G(1,:); B(2,3:3:end)=G(2,:);
        B(3,2:3:end)=G(2,:); B(3,3:3:end)=G(1,:);
        S=J\[(1-eta)/2*bottom+(1+eta)/2*top; ...
            (1-xi)/2*left+(1+xi)/2*right];
        H=zeros(3,12);
        for k=1:3, H(k,k:3:end)=N; end
        K=K+(B'*Db*B+S'*Ds*S)*determinant;
        M=M+rho*H'*diag([t,t^3/12,t^3/12])*H*determinant;
        P(1:3:end)=P(1:3:end)+N'*determinant;
        area=area+determinant;
    end
end
K=(K+K')/2; M=(M+M')/2;
end

function row=tied(xi,eta,direction,X)
[N,dN]=plate.shape(xi,eta); tangent=dN(direction,:)*X;
row=zeros(1,12); row(1:3:end)=dN(direction,:);
row(2:3:end)=-tangent(1)*N; row(3:3:end)=-tangent(2)*N;
end
