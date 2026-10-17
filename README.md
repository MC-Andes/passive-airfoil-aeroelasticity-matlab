# Passive Airfoil Aeroelasticity MATLAB

[![Verify](https://github.com/MC-Andes/passive-airfoil-aeroelasticity-matlab/actions/workflows/verify.yml/badge.svg)](https://github.com/MC-Andes/passive-airfoil-aeroelasticity-matlab/actions/workflows/verify.yml)

Solver propio MATLAB para el álabe original **fijo en una corriente de aire**:
FEM sólido Tet4/Tet10, transferencia conservativa, VLM/UVLM, acoplamiento
estático y dinámico, y una rama Navier–Stokes laminar 2D/3D para estela viscosa.
Gmsh/OpenCASCADE se usan únicamente para preparar geometría y mallas.

El estudio del CAD conserva el sólido y usa montaje/material/slider como
[hipótesis explícitas](cases/cad_original/assumptions.json), autorizadas ante la
ausencia de datos experimentales. La nueva familia prismática converge en
observables TE y la banda estructural declarada; las dos mallas históricas
siguen siendo un diagnóstico de falta de convergencia. Ver el
[informe de resultados y alcance](docs/cad-original-results.md).

El [portal MC-Andes](https://mc-andes.github.io/projects/) reúne los proyectos
de la organización. La [guía de publicación](docs/portal-publication.md) registra
el alcance y la reproducción de los recursos de esta ficha.

## Inicio rápido

Requiere **MATLAB R2026a**; el núcleo usa MATLAB base. No hacen falta PDE Toolbox,
Python ni un solver externo para ejecutar las mallas incluidas.

```bash
git clone https://github.com/MC-Andes/passive-airfoil-aeroelasticity-matlab.git
cd passive-airfoil-aeroelasticity-matlab
```

En MATLAB, desde la raíz del proyecto:

```matlab
addpath('matlab')
report = run_verification();       % siete suites: estructura, CAD, UVLM/FSI y Navier–Stokes
run('matlab/run_fem_demo.m')        % ejemplo paramétrico, modos y vibración libre
comparison = run_tail_benchmarks(); % nueve geometrías PLA publicadas
result = run_tail_response();       % deformación TE y presión oscilante impuesta
```

Desde una terminal donde `matlab` esté en `PATH`:

```bash
matlab -batch "addpath('matlab'); run_verification(); run('matlab/run_fem_demo.m');"
```

Se generan CSV, figuras y archivos MAT en `results/`, excluido de Git. Los
scripts resuelven las rutas desde el proyecto y no necesitan archivos personales.

La ruta nueva de [borde de salida](docs/trailing-edge.md) incluye comparación
con [casos publicados](docs/case-catalog.md). Mantiene el modelo desacoplado:
la oscilación de presión del ejemplo es impuesta. La [ruta CAD acoplada](docs/cad-coupled.md) calcula realimentación y estela
circulatoria; la [rama viscosa](docs/viscous-flow.md) calcula presión, déficit
y vorticidad bajo sus propios supuestos y Reynolds.

## Ejemplos disponibles

| Entrada | Qué ejecuta |
| --- | --- |
| `matlab/run_fem_demo.m` | Sólido NACA paramétrico, cargas VLM de una vía, modos y Newmark estructural. |
| `matlab/run_fem_reference.m` | FEM sobre malla derivada del CAD, con material y empotramiento supuestos. |
| `matlab/run_fem_convergence.m` | Diagnóstico de dos mallas; informa falta de convergencia. |
| `matlab/inspect_reference.m` | Inspección de triangulación y envolventes aerodinámicas. |
| `matlab/run_demo.m` | Modelo reducido ilustrativo con barrido y transitorio cuasiestacionario. |
| `matlab/run_reference.m` | Superficie derivada del CAD con estructura reducida ilustrativa. |
| `matlab/run_verification.m` | Siete suites, reportes y resumen de verificación. |
| `matlab/run_cad_original.m` | Campaña estructural del CAD con cargas/sensores físicos y criterios explícitos. |
| `matlab/run_cad_aeroelastic.m` | Equilibrios de una/dos vías y perturbación libre con UVLM del CAD. |
| `matlab/run_cad_mechanisms.m` | Hipótesis A/B/C de guía sobre el sólido conectado. |
| `matlab/run_cad_viscous.m` | Estela laminar 2D/3D y pared móvil prescrita en el exterior CAD. |
| `matlab/run_tail_benchmarks.m` | Nueve geometrías PLA publicadas: frecuencias adimensionales, errores y refinación. |
| `matlab/run_tail_response.m` | Placa equivalente: cola estática, presión temporal impuesta y respuesta en frecuencia. |

Para ejecutar el caso de referencia:

```matlab
run('matlab/inspect_reference.m')
run('matlab/run_fem_reference.m')
run('matlab/run_fem_convergence.m')
```

La malla refinada requiere más tiempo y memoria que el ejemplo paramétrico.
Los scripts sobrescriben las salidas homónimas al repetirse.

## Qué está implementado

- Placa MITC4 propia para flexión delgada, masa consistente y presión distribuida.
- Reducción modal desde matrices FEM de placas o sólidos; presión temporal y observables comunes.
- Comparación reproducible de nueve geometrías PLA con datos medidos atribuidos (27 observaciones).
- Elementos sólidos Tet4/Tet10 con rigidez y masa consistente calculadas desde cero.
- Ensamblaje disperso, materiales isotrópicos u ortótropos y regiones materiales.
- Eliminación de desplazamientos prescritos, reacciones, tensiones y energía.
- Modos estructurales normalizados en masa y Newmark de aceleración media.
- VLM de envergadura finita con fuerzas de Kutta–Joukowski.
- Interpolación afín y transferencia por trabajo virtual conservando fuerza y momento.
- Datos de referencia de Tough PLA y ABS impresos, con procedencia y supuestos explícitos.
- Modelo reducido independiente para estudios cuasiestacionarios y callbacks propios.

El FEM calcula sus matrices desde la malla y los materiales; no delega el
cálculo a Gmsh o a PDE Toolbox. Gmsh/OpenCASCADE solo se utiliza en las
herramientas opcionales para generar entradas geométricas.

## Hipótesis y límites

La ruta principal usa elasticidad pequeña con Tet10; existe una formulación
no lineal total Lagrangiana verificada para comparar cuando la deformación lo
requiera. Las guías, contactos afines y Coulomb son hipótesis de pequeñas
excursiones. El STEP sigue siendo un sólido conectado; no se inventan bisagras.
No se identifica automáticamente el montaje ni el mecanismo fabricado.

VLM/UVLM representan flujo adherido inviscido. UVLM incorpora memoria de
circulación y fuerzas no estacionarias; `pressureJump` sigue siendo una carga
equivalente delgada, no Cp de ambas pieles. La rama Navier–Stokes resuelve casos
laminares declarados y necesita su propio estudio de resolución/penalización.
No incorpora modelos calibrados de transición o turbulencia. **No hay una
velocidad de flutter validada experimentalmente** ni datos PIV del espécimen.

La ley de referencia adopta Tough PLA con E=2,617 GPa y rho=1160 kg/m³ o ABS
con E=1,813 GPa y rho=1084 kg/m³, a partir de
[Malim et al. (2023)](https://doi.org/10.3390/aerospace10030256). En ambos casos
nu=0,35 es un supuesto. Esas propiedades no caracterizan automáticamente la
pieza fabricada ni cualquier orientación o porcentaje de relleno.

## Documentación completa

Comenzar por el [índice de documentación](docs/index.md).

| Tema | Documento |
| --- | --- |
| Instalación, comandos y solución de problemas | [Inicio rápido](docs/getting-started.md) |
| Arquitectura, ejes, unidades y grados de libertad | [Arquitectura](docs/architecture.md) |
| Formulación estructural y transferencia | [FEM propio](docs/solid-fem.md) |
| Funciones, argumentos y salidas | [Referencia de API](docs/api-reference.md) |
| Ejemplos e interpretación de resultados | [Ejemplos](docs/examples.md) |
| Formatos, procedencia y preparación de otro CAD | [Datos](docs/data-and-preprocessing.md) |
| Pruebas, convergencia y validación | [Verificación](docs/verification.md) |
| Modelo reducido y sus extensiones | [Modelo reducido](docs/reduced-model.md) |
| Información necesaria del laboratorio | [Datos experimentales](docs/experimental-data.md) |
| Bibliografía y materiales | [Referencias](docs/references.md) |
| Acoplamiento sobre el CAD | [Mecánica y UVLM/FSI](docs/cad-coupled.md) |
| Déficit, presión y vorticidad viscosos | [Navier–Stokes](docs/viscous-flow.md) |
| Mejoras pendientes | [Hoja de ruta](docs/roadmap.md) |

## Contenido publicado

Se incluyen únicamente código, pruebas, configuración de verificación,
documentación y las entradas numéricas necesarias para correr los ejemplos.
Las mallas inicial/refinada, triangulación, secciones y metadatos están en
`data/`; sus hashes están en [el manifiesto](data/manifest.json).

El CAD original, STL redundantes, miniaturas, artículos completos, resultados,
dependencias descargadas y entornos locales quedan fuera del repositorio. El
CAD original no hace falta para ejecutar los CSV preparados. Las herramientas
opcionales aceptan un STEP propio si se desea generar otra geometría.

Para comprobar integridad y documentación sin MATLAB:

```bash
python3 tools/check_repository.py
```

## Desarrollo y licencia

Ver [CONTRIBUTING](CONTRIBUTING.md), [CHANGELOG](CHANGELOG.md) y
[CITATION.cff](CITATION.cff). El código y la documentación se publican con
[licencia MIT](LICENSE), siguiendo la licencia por defecto de MC-Andes. Las
dependencias externas conservan sus licencias. El extracto experimental
Huang–Zhang se distribuye bajo CC BY 4.0, con atribución y transformaciones en
[avisos de terceros](THIRD_PARTY_NOTICES.md).

La importación inicial establece `main` como excepción de arranque a la regla
de integración por PR de [MC-Andes](https://github.com/MC-Andes/standards).
Los cambios posteriores se revisan mediante PR; no se declara una aprobación
independiente ni una validación física para esta importación inicial.

## Desarrollo sobre el CAD original

El [plan CAD](docs/plan-cad-original.md) tiene una
[campaña estructural ejecutable](docs/cad-original.md) y
[contratos de caso](cases/cad_original/case.json). Incorpora cargas de presión
sobre un parche físico, nueve sensores materiales cerca del borde de salida,
atlas geométrico, Tet10/Tri6, escenarios A/B/C, UVLM/FSI y flujo viscoso
laminar propio. `run_verification()` reúne ahora 69 controles en siete suites.
El [informe](docs/cad-original-results.md) cuantifica convergencia y límites.
Los resultados siguen siendo condicionales a material, apoyos y guía supuestos. El objetivo principal y el núcleo MATLAB propio
se mantienen.
