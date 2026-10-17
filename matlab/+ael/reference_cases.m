function cases = reference_cases()
%REFERENCE_CASES Geometry-keyed published plate cases; material is a nominal model.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
T=readtable(fullfile(root,'data','benchmarks','huang2026_modal.csv'),'TextType','string');
P=jsondecode(fileread(fullfile(root,'data','benchmarks','huang2026_provenance.json')));
ids=unique(T.case_id,'stable'); cases=struct([]);
for k=1:numel(ids)
    rows=T(T.case_id==ids(k),:); c=struct(); c.id=ids(k);
    c.chord=P.nominal_geometry.chord_m; c.span=rows.aspect_ratio(1)*c.chord;
    c.flexibleFraction=rows.flexible_chord_fraction(1);
    c.rigidThickness=P.nominal_geometry.rigid_thickness_m;
    c.flexibleThickness=P.nominal_geometry.flexible_thickness_m;
    c.material=struct('E',P.nominal_model_parameters.E_Pa, ...
        'nu',P.nominal_model_parameters.nu,'rho',P.nominal_model_parameters.rho_kg_m3);
    c.reference=rows; c.sourceDOI=P.dataset_doi;
    c.scope='nominal flat plate, experimental dimensionless modal comparison';
    if k==1, cases=c; else, cases(k)=c; end %#ok<AGROW>
end
end
