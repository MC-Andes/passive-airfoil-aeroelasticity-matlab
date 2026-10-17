function [point,N,r,ok] = face_point(mesh,face,xy,seed)
%FACE_POINT Invert a Tri3/Tri6 graph face in physical x/span coordinates.
X=mesh.nodes(face,:);nc=numel(face);tol=1e-11*max(1,max(abs(X),[],'all'));
J=(X(2:3,1:2)-X(1,1:2))';
if nargin<4
    if rcond(J)>1e-12,r=(J\(xy(:)-X(1,1:2)'))';else,r=[1/3,1/3];end
else,r=seed(:)';end
for iteration=1:30
    [N,G]=fem3d.tri_shape(r,nc);point=N*X;defect=point(1:2)-xy(:)';
    if norm(defect,inf)<tol,break;end
    tangent=X(:,1:2)'*G;if rcond(tangent)<1e-13,break;end
    delta=(tangent\defect')';accepted=false;
    for search=0:10
        candidate=r-2^(-search)*delta;sh=fem3d.tri_shape(candidate,nc);
        if norm(sh*X(:,1:2)-xy(:)',inf)<norm(defect,inf),r=candidate;accepted=true;break;end
    end
    if ~accepted,break;end
end
[N,~]=fem3d.tri_shape(r,nc);point=N*X;L=[1-sum(r),r];
ok=norm(point(1:2)-xy(:)',inf)<tol && all(L>=-1e-8) && all(L<=1+1e-8);
end
