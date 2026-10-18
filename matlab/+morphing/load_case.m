function c = load_case(root,path)
%LOAD_CASE Updated PETG/servo study; distinct from historical passive cases.
if nargin<2,path=fullfile(root,'cases/morphing_wing/case.json');end
c=jsondecode(fileread(path));assert(c.schemaVersion==2,'morphing:Case','Require updated study schema 2.');
a=jsondecode(fileread(fullfile(root,c.geometry.audit)));e=jsondecode(fileread(fullfile(root,c.studyData)));
assert(strcmp(a.source_step_sha256,c.geometry.sourceSha256) && strcmp(e.geometry.updatedStepSha256,c.geometry.sourceSha256) && a.step.brep_valid,'morphing:Hash','Geometry/source contract mismatch.');
assert(isequal(c.geometry.axisOrder(:)',[1,3,2]) && max(abs(c.geometry.origin_m(:)-a.frame.origin_m(:)))<1e-12,'morphing:Frame','Inconsistent physical coordinates.');
assert(abs(c.geometry.referenceArea_m2-c.geometry.chord_m*c.geometry.span_m)<1e-12,'morphing:Reference','Half-wing reference area mismatch.');
assert(all(c.flow.speeds_m_s>0) && c.flow.rho_kg_m3>0 && c.flow.mu_Pa_s>0 && c.slider.positions==11,'morphing:Case','Invalid flow/slider data.');
assert(c.material.E_Pa>0 && c.material.rho_kg_m3>0 && c.material.nu>0 && c.material.nu<.5,'morphing:Material','Invalid conditional material.');
assert(c.supports.insertion_m>0 && c.supports.insertion_m<=c.geometry.span_m,'morphing:Support','Invalid shaft insertion.');
assert(isequal(sort(c.supports.cadFaces(:)),sort(a.shaft_socket.faces(:))) && isequal(sort(c.slider.cadFaces(:)),sort(a.slider_drive.faces(:))),'morphing:Labels','Source-bound interfaces differ from audit.');
assert(abs(c.geometry.chord_m-e.geometry.chord_m)<1e-12 && abs(c.geometry.span_m-e.geometry.span_m)<1e-12 && abs(c.flow.rho_kg_m3-e.flow.rho_kg_m3)<1e-12 && abs(c.flow.mu_Pa_s-e.flow.mu_Pa_s)<1e-12 && isequal(c.flow.speeds_m_s,e.flow.speeds_m_s) && isequal(c.flow.alpha_deg,e.flow.alpha_deg),'morphing:Study','Study data and operational flow/geometry disagree; define an explicit variant instead.');
c.audit=a;c.study=e;c.caseSha256=cad.sha256(path);
end
