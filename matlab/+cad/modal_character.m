function report = modal_character(mesh,modes,config)
%MODAL_CHARACTER Consistent modal kinetic energy by region and direction.
% Posterior region uses whole elements whose material centroid has x/c>=.70;
% it is a declared regional observable, not a stress or joint classifier.
T=mesh.tetra;X=mesh.nodes;Phi=modes.Phi;nm=size(Phi,2);mass=zeros(nm,3);tail=mass;
reference=fem3d.elevate(fem3d.box_mesh([1,1,1],[1,1,1]));xe=reference.nodes(reference.tetra(1,:),:);[~,me,~,vol]=fem3d.tet10(xe,eye(6),1);massReference=me(1:3:end,1:3:end)/vol;
for start=1:2500:size(T,1)
 ix=start:min(start+2499,size(T,1));tt=T(ix,:);ne=numel(ix);a=X(tt(:,2),:)-X(tt(:,1),:);b=X(tt(:,3),:)-X(tt(:,1),:);c=X(tt(:,4),:)-X(tt(:,1),:);volume=abs(sum(a.*cross(b,c,2),2))/6;center=reshape(mean(reshape(X(tt(:,1:4)',1),4,ne),1),[],1);posterior=center>=.7*config.geometry.chord_m;
 for axis=1:3
  U=reshape(Phi(3*tt'-3+axis,:),10,ne,nm);U=permute(U,[1,3,2]);energy=squeeze(sum(U.*pagemtimes(massReference,U),1))';
  mass(:,axis)=mass(:,axis)+sum(volume.*energy,1)';tail(:,axis)=tail(:,axis)+sum((volume.*posterior).*energy,1)';
 end
end
fraction=mass./sum(mass,2);posterior=sum(tail,2)./sum(mass,2);frequency=modes.frequencyHz;
report=table((1:nm)',frequency,posterior,fraction(:,1),fraction(:,2),fraction(:,3),'VariableNames',{'mode','frequency_Hz','posterior_modal_mass_fraction','chord_fraction','span_fraction','height_fraction'});
end
