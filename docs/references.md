# Referencias y uso dentro del proyecto

La revisión es enfocada y no pretende ser una búsqueda sistemática exhaustiva.
Las propiedades publicadas se usan como referencia inicial, no como una
caracterización de la pieza del usuario. No se redistribuyen los artículos.

## Morphing pasivo y validación aerodinámica

- **Taguchi et al. (2020).** Experimental Study about the Deformation and
  Aerodynamic Characteristics of the Passive Morphing Airfoil. Transactions of
  the Japan Society for Aeronautical and Space Sciences, 63(1), 18–23.
  [DOI](https://doi.org/10.2322/tjsass.63.18).
  Motiva comparar forma y fuerzas, y distinguir zona anterior rígida y zona flexible.
- **Kai et al. (2023; online 2022).** Low-Speed Wind Tunnel Testing of a Passive
  Camber Morphing Airfoil Using a 3D-Printed Compliant Mechanism. APISAT 2021,
  LNEE 912, 1019–1030.
  [DOI](https://doi.org/10.1007/978-981-19-2689-1_78).
  Referencia de arquitectura compliant impresa y pruebas de NACA0012 pasivo.
- **Ochi et al. (2024).** Aeroelastic Simulation and Experimental Validation
  of a 3D-Printed Passive Morphing Airfoil. AIAA Journal.
  [DOI](https://doi.org/10.2514/1.J063000).
  El resumen consultado describe comparación Euler/RANS y validación. No se
  extrajeron de este trabajo propiedades numéricas para los materiales del código.

## Manufactura aditiva y estructura

- **Malim, Mourousias, Marinus y De Troyer (2023).** Structural Design of a
  Large-Scale 3D-Printed High-Altitude Propeller: Methodology and Experimental
  Validation. Aerospace, 10(3), 256.
  [DOI](https://doi.org/10.3390/aerospace10030256).
  Su tabla 3 es la fuente de E, densidad y resistencias últimas almacenadas para
  Tough PLA/ABS. Poisson 0,35 y el amortiguamiento no se atribuyen a esa tabla.
- **Dogea, Yan y Millar (2023).** Testing and evaluation of the structural
  performance of a 3D-printed polylactic acid aircraft wing rib.
  Discover Mechanical Engineering, 2, 6.
  [DOI](https://doi.org/10.1007/s44245-023-00013-8).
  Referencia de sensibilidad de una nervadura impresa a carga y montaje.
- **Gonabadi, Yadav y Bull (2020).** The effect of processing parameters on the
  mechanical characteristics of PLA produced by a 3D FFF printer.
  The International Journal of Advanced Manufacturing Technology, 111, 695–709.
  [DOI](https://doi.org/10.1007/s00170-020-06138-4).
  Fundamenta distinguir orientación y propiedades efectivas de impresión.
- **Huang y Zhang (2026).** Comparative modal analysis of beam-plate coupled
  models for isotropic active camber morphing wings.
  Aerospace Science and Technology, 178, 112449.
  [DOI](https://doi.org/10.1016/j.ast.2026.112449).
  Orienta la selección de modelos para flexión en cuerda/envergadura y torsión;
  no se afirma haber implementado su modelo mixto.
- **Huang y Zhang, dataset versión 2 (2026).** Experimental and numerical data
  from comparative analysis of candidate structural models for active camber
  morphing wings. [Repositorio](https://data.mendeley.com/datasets/8jxw5533py/2).
  Referencia candidata para validación posterior. Sus ensayos no se han
  incorporado al solver ni ajustado a la geometría de referencia.

## Métodos y herramientas

- [NASA, disposición de mallas VLM](https://ntrs.nasa.gov/citations/19760021094):
  contexto metodológico para el método de red de vórtices.
- [SHARPy](https://doi.org/10.21105/joss.01885): referencia de arquitectura
  aeroelástica; este proyecto no reproduce sus capacidades UVLM.
- [Bathe, placas y láminas](https://ocw.mit.edu/courses/res-2-002-finite-element-procedures-for-solids-and-structures-spring-2010/resources/lecture-20/):
  contexto para formulaciones de lámina y bloqueo por cortante.
- [Gmsh](https://gmsh.info/doc/texinfo/): preparación geométrica opcional;
  no calcula las matrices o soluciones estructurales de este repositorio.
- [OpenCASCADE BRepGProp](https://dev.opencascade.org/doc/refman/html/class_b_rep_g_prop.html):
  integración adaptativa y no adaptativa de propiedades geométricas.
- [MathWorks MATLAB Actions](https://github.com/matlab-actions/setup-matlab):
  instalación de MATLAB en CI pública. La licencia del código no concede
  una licencia de MATLAB para una instalación local.

La formulación implementada y sus convenciones se describen en
[FEM sólido](solid-fem.md), [arquitectura](architecture.md) y [API](api-reference.md).
Al publicar resultados científicos derivados, registrar también el commit
del código y las referencias de los modelos y materiales realmente utilizados.
