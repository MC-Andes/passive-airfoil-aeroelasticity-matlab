function hit = ray(mesh,xy,faceIds)
%RAY Vertical intersections with a faceted solid boundary in canonical axes.
% Tri6 intersections use exact inverse isoparametric surface mapping.
if nargin<3, faceIds=(1:size(mesh.faces,1))'; end
if size(mesh.faces,2)==6
    hit=curved_ray(mesh,xy,faceIds);return;
end
faces=mesh.faces(faceIds,:); [tri,map]=fem3d.plot_faces(faces); X=mesh.nodes;
a=X(tri(:,1),:); b=X(tri(:,2),:); c=X(tri(:,3),:);
tol=1e-11*max(1,max(abs(X),[],'all'));
take=xy(1)>=min([a(:,1),b(:,1),c(:,1)],[],2)-tol & ...
    xy(1)<=max([a(:,1),b(:,1),c(:,1)],[],2)+tol & ...
    xy(2)>=min([a(:,2),b(:,2),c(:,2)],[],2)-tol & ...
    xy(2)<=max([a(:,2),b(:,2),c(:,2)],[],2)+tol;
idx=find(take); a=a(idx,:); b=b(idx,:); c=c(idx,:);
e=b(:,1:2)-a(:,1:2); g=c(:,1:2)-a(:,1:2); d=xy-a(:,1:2); det=e(:,1).*g(:,2)-e(:,2).*g(:,1);
v=(d(:,1).*g(:,2)-d(:,2).*g(:,1))./det;
w=(e(:,1).*d(:,2)-e(:,2).*d(:,1))./det;
inside=abs(det)>tol^2 & v>=-1e-9 & w>=-1e-9 & v+w<=1+1e-9;
idx=idx(inside); bary=[1-v(inside)-w(inside),v(inside),w(inside)];
z=sum(bary.*[a(inside,3),b(inside,3),c(inside,3)],2);
[z,sortIndex]=sort(z); idx=idx(sortIndex); bary=bary(sortIndex,:);
keep=[true;diff(z)>tol]; if isempty(z), keep=false(0,1); end
hit=struct('z',z(keep),'face',faceIds(map(idx(keep))),'triangle',tri(idx(keep),:),'barycentric',bary(keep,:),'parametric',bary(keep,2:3));
end

function hit=curved_ray(mesh,xy,faceIds)
F=mesh.faces(faceIds,:);X=mesh.nodes;tol=1e-11;
control=zeros(size(F,1),6,2);
for k=1:3,control(:,k,:)=reshape(X(F(:,k),1:2),[],1,2);end
edges=[1,2;2,3;3,1];
for k=1:3,control(:,k+3,:)=reshape(2*X(F(:,k+3),1:2)-(X(F(:,edges(k,1)),1:2)+X(F(:,edges(k,2)),1:2))/2,[],1,2);end
lo=reshape(min(control,[],2),[],2);hi=reshape(max(control,[],2),[],2);
use=find(all(xy>=lo-tol & xy<=hi+tol,2) & all(hi-lo>tol,2));
z=zeros(0,1);ids=zeros(0,1);r=zeros(0,2);
for k=use'
    [p,~,param,ok]=cad.face_point(mesh,F(k,:),xy);
    if ~ok
        for seed=[.1,.1;.8,.1;.1,.8]'
            [p,~,param,ok]=cad.face_point(mesh,F(k,:),xy,seed);
            if ok,break;end
        end
    end
    if ok,z(end+1,1)=p(3);ids(end+1,1)=faceIds(k);r(end+1,:)=param;end %#ok<AGROW>
end
[z,order]=sort(z);ids=ids(order);r=r(order,:);keep=[true;diff(z)>tol];if isempty(z),keep=false(0,1);end
hit=struct('z',z(keep),'face',ids(keep),'parametric',r(keep,:), ...
    'barycentric',[1-sum(r(keep,:),2),r(keep,:)],'triangle',mesh.faces(ids(keep),1:3));
end
