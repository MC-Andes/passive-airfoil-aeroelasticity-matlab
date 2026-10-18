function estimate = resource_viscous(config)
%RESOURCE_VISCOUS Transparent uniform-grid budget, not a DNS capability claim.
Re=config.flow.rho_kg_m3*config.flow.speed_m_s*config.geometry.chord_m/config.flow.mu_Pa_s;
delta=5/sqrt(Re);spacing=delta/10;box=[7,3.5,3];count=ceil(box/spacing);cells=prod(count);
estimate=struct('Re',Re,'laminarFlatPlateDelta99_over_c',delta,'cellsAcrossDelta',10,'uniformCellWidth_over_c',spacing,'domain_over_c',box,'gridCount',count,'cells',cells,'minimumEightScalarArrays_bytes',8*8*cells,'workingThirtyTwoArrays_bytes',32*8*cells,'hostRAM_bytes',24*1024^3,'fitsHost',32*8*cells<24*1024^3,'scope','order-of-magnitude uniform-grid laminar boundary-layer budget; curved immersed wall, turbulence scales and pressure iterations can require much more','decision','retain verified low-Re laminar studies and attached-flow UVLM; do not claim tunnel-Re DNS or calibrated URANS/LES');
end
