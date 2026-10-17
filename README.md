# Passive Airfoil Aeroelasticity MATLAB

[![Verify](https://github.com/MC-Andes/passive-airfoil-aeroelasticity-matlab/actions/workflows/verify.yml/badge.svg)](https://github.com/MC-Andes/passive-airfoil-aeroelasticity-matlab/actions/workflows/verify.yml)

Solver de desarrollo para perfiles 3D **fijos en una corriente de aire** que se
deforman pasivamente. Implementa un FEM sólido propio en MATLAB, un VLM de
herraduras y transferencia conservativa de cargas. La ruta principal es
**desacoplada de una vía**: aerodinámica sobre la forma inicial → cargas → FEM.

**Estado:** prototipo de investigación, con verificación numérica y sin
validación experimental de la pieza de referencia. Su respuesta estructural
todavía no converge entre las dos mallas suministradas: el desplazamiento
máximo cambia aproximadamente **66,1%**. Los resultados de ese caso no son
predicciones cuantitativas validadas del álabe ni de su slider.

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
report = run_verification();       % 33 controles: VLM/reducido, sólido y cola
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
la oscilación de presión del ejemplo es impuesta. La estela no estacionaria y
el acoplamiento en dos vías quedan como etapas siguientes.

## Ejemplos disponibles

| Entrada | Qué ejecuta |
| --- | --- |
| `matlab/run_fem_demo.m` | Sólido NACA paramétrico, cargas VLM de una vía, modos y Newmark estructural. |
| `matlab/run_fem_reference.m` | FEM sobre malla derivada del CAD, con material y empotramiento supuestos. |
| `matlab/run_fem_convergence.m` | Diagnóstico de dos mallas; informa falta de convergencia. |
| `matlab/inspect_reference.m` | Inspección de triangulación y envolventes aerodinámicas. |
| `matlab/run_demo.m` | Modelo reducido ilustrativo con barrido y transitorio cuasiestacionario. |
| `matlab/run_reference.m` | Superficie derivada del CAD con estructura reducida ilustrativa. |
| `matlab/run_verification.m` | Tres suites, reportes y resumen de verificación. |
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
- Elementos sólidos Tet4 con rigidez y masa consistente calculadas desde cero.
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

El FEM sólido es lineal de pequeña deformación y rotación. Tet4 puede resultar
demasiado rígido en flexión de pieles delgadas y en materiales casi
incompresibles. No incluye contacto, holgura, fricción, hiperelasticidad ni la
cinemática real del slider. El ejemplo del CAD considera un sólido unido y
empotra el 30% anterior de cuerda en toda la envergadura; son hipótesis de desarrollo.

El VLM es inviscido, de superficie delgada y con estela recta cuasiestacionaria.
No resuelve transición, separación, resistencia viscosa, memoria de estela ni
masa añadida. `pressureJump` no es Cp real de las dos caras y `CDi` no es
resistencia total. **No hay una predicción de flutter validada**. Los modos
del FEM y los autovalores del modelo reducido tienen alcances distintos.

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
