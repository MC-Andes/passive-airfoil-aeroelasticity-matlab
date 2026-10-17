function comparisons = run_cad_modal_shape_comparison(campaign)
%RUN_CAD_MODAL_SHAPE_COMPARISON MAC over chord/span, not TE signatures alone.
root=fileparts(fileparts(mfilename('fullpath')));previous=[];comparisons=cell(2,1);stations=[.4,.55,.7,.85,.999];
for level=1:3
 s=load(fullfile(campaign,"tet10_level_"+level,'solution.mat'),'config','modes');c=s.config;mesh=fem3d.elevate(fem3d.read_mesh(fullfile(root,string(c.discretization.meshes{level})),c.geometry.axisOrder,c.geometry.origin_m));H=[];
 for x=stations,c.sensors.x_over_c=x;sensor=cad.sensors(mesh,c);H=[H;sensor.H];end %#ok<AGROW>
 observation=kron(H,speye(3))*s.modes.Phi;eta=c.sensors.eta(:);wy=diff([0;(eta(1:end-1)+eta(2:end))/2;1]);wx=diff([stations(1);((stations(1:end-1)+stations(2:end))/2)';stations(end)]);weights=repelem(kron(wx(:),wy),3);current=struct('observed',observation,'frequencyHz',s.modes.frequencyHz);
 if level>1
  A=previous.observed;B=observation;MAC=abs(A'*(weights.*B)).^2./(sum(weights.*abs(A).^2,1)'*sum(weights.*abs(B).^2,1));comparisons{level-1}=struct('from',level-1,'to',level,'MACmatrix',MAC,'pairedMAC',diag(MAC),'stations_x_over_c',stations,'eta',eta,'weightDefinition','tensor trapezoid chord/span weights on xyz; exposed upper material surface','frequencyRelativeChange',abs(current.frequencyHz./previous.frequencyHz-1));
 end
 previous=current;
end
output=fullfile(campaign,'modal_shape_comparison');if ~isfolder(output),mkdir(output);end;cad.write_json(fullfile(output,'comparison.json'),comparisons);save(fullfile(output,'comparison.mat'),'comparisons');disp(comparisons{end});
end
