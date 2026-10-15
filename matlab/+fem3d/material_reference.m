function material = material_reference(name)
%MATERIAL_REFERENCE Published AM properties, NOT properties of the supplied CAD.
% Poisson .35 for Malim datasets is an explicit assumption; not measured there.
arguments
    name (1,1) string = "tough_pla_malim2023"
end
material=struct('name',name,'nu',.35,'rho',NaN,'E',NaN);
switch name
    case "tough_pla_malim2023"
        material.E=2617e6; material.rho=1160;
        material.ultimateStrengthXYZ=[47,47,38]*1e6;
        material.source="Malim et al. (2023), Table 3; DOI 10.3390/aerospace10030256";
        material.url="https://doi.org/10.3390/aerospace10030256";
        material.printing="FDM, 100% infill, raster -45/+45 deg, layer 0.2 mm";
        material.assumptions="nu=0.35 assumed; isotropic small-strain homogenization; no strength criterion calibrated";
    case "abs_malim2023"
        material.E=1813e6; material.rho=1084;
        material.ultimateStrengthXYZ=[37.5,37.5,20]*1e6;
        material.source="Malim et al. (2023), Table 3; DOI 10.3390/aerospace10030256";
        material.url="https://doi.org/10.3390/aerospace10030256";
        material.printing="FDM, 100% infill, raster -45/+45 deg, layer 0.2 mm";
        material.assumptions="nu=0.35 assumed; isotropic small-strain homogenization; no strength criterion calibrated";
    otherwise
        error('fem3d:Material','Unknown reference material: %s',name);
end
material.D=fem3d.isotropic(material.E,material.nu);
end
