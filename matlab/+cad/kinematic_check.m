function report = kinematic_check(mesh,u,config)
%KINEMATIC_CHECK Geometric-linearity review using local displacement gradients.
T=mesh.tetra;X=mesh.nodes;U=reshape(u,3,[])';rotation=zeros(size(T,1),1);strain=rotation;
if size(T,2)==10,[~,G]=fem3d.tet10_shape([.25,.25,.25]);else,G=[-1,-1,-1;eye(3)];end
for e=1:size(T,1)
    xe=X(T(e,:),:);xe=xe-xe(1,:);ue=U(T(e,:),:);ue=ue-ue(1,:);grad=G/(xe'*G);H=ue'*grad;
    rotation(e)=norm([H(3,2)-H(2,3),H(1,3)-H(3,1),H(2,1)-H(1,2)])/2;
    strain(e)=norm((H+H')/2,'fro');
end
flexibleLength=(1-config.supports.x_over_c)*config.geometry.chord_m;
relative=max(vecnorm(U,2,2))/flexibleLength;
report=struct('maximumDisplacementOverFlexibleLength',relative,'maximumSmallRotation_rad',max(rotation), ...
    'maximumCentroidInfinitesimalStrainNorm',max(strain),'comparisonTriggerDisplacement',.05, ...
    'comparisonTriggerRotation_rad',.1,'nonlinearComparisonTriggered',relative>=.05 || max(rotation)>=.1, ...
    'interpretation','review triggers, not universal validity bounds; local singular stress is not calibrated strength');
end
