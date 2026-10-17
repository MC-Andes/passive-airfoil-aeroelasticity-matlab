# Verificación, convergencia y validación

## Ejecutar los controles

```matlab
addpath('matlab')
report = run_verification();
```

El wrapper ejecuta `run_tests`, `run_fem_tests` y `run_tail_tests`. Las suites usan referencias
analíticas o invariantes de física y fallan con un error cuando no cumplen sus
criterios. Los resultados se guardan bajo `results/verification/` y el FEM
también mantiene `results/fem_test_results.csv` por compatibilidad.

## FEM sólido: doce controles

| Control | Referencia / propiedad contrastada |
| --- | --- |
| `tetra_rigid_body` | Seis modos rígidos y deformación nula bajo traslación/rotación infinitesimal. |
| `consistent_mass` | Energía cinética de traslación igual a rho·V·velocidad². |
| `orientation_invariance` | Invertir el orden de nodos no cambia las matrices físicas ensambladas. |
| `affine_patch_test` | Campo afín exacto con tracciones derivadas de su tensor de tensiones y apoyos no nulos. |
| `pressure_balance` | Presión uniforme sobre superficie cerrada sin fuerza ni momento resultantes. |
| `cantilever_refinement` | Tendencia de flexión frente a Euler–Bernoulli en tres refinamientos. |
| `modal_residual` | Residual del problema propio y ortogonalidad de masa. |
| `newmark_free_mode` | Historia cosenoidal exacta del primer modo y conservación de energía sin amortiguamiento. |
| `conservative_transfer` | Reproducción afín, fuerza, momento y trabajo virtual. |
| `material_orthotropic_limit` | Recuperación de la matriz isotrópica al igualar constantes direccionales. |
| `unsupported_mesh_rejected` | Rechazo de un sólido con movimientos rígidos sin restringir. |
| `cad_volume` | Integral de volumen de los tetraedros frente a auditoría independiente y compatibilidad de superficie. |

El test de voladizo acepta menos de 15% de error en su última malla y exige
reducción monótona del error. Es un control de tendencia y de la formulación,
no una promesa de precisión de 15% para otras geometrías. En la referencia
local los errores fueron aproximadamente 70,4%, 39,3% y 14,9%. Esto expone la
rigidez artificial de Tet4 en flexión cuando el refinamiento es insuficiente.

## VLM y modelo reducido: once controles

| Control | Referencia / propiedad contrastada |
| --- | --- |
| `rigid_symmetry` | Sustentación antisimétrica con alpha, carga cero a incidencia cero y residual de no penetración. |
| `thin_airfoil_limit` | Pendiente próxima a 2π para gran alargamiento y ángulo pequeño. |
| `coefficient_scaling` | Escala dimensional rho·U² y coeficientes adimensionales consistentes. |
| `mesh_convergence` | Diferencia de sustentación entre dos discretizaciones de un ala de referencia. |
| `virtual_work` | Fuerzas generalizadas frente a trabajo de una cinemática rotacional. |
| `analytical_feedback` | Equilibrio y linealización de un modelo de realimentación soluble. |
| `divergence_failure` | Rechazo de equilibrio singular en el modelo analítico de prueba. |
| `free_oscillator` | Oscilador amortiguado exacto y decaimiento de energía. |
| `steady_transient` | Conservación de un equilibrio en la integración temporal. |
| `affine_modal_interpolation` | Reproducción de un modo afín interpolado. |
| `reference_geometry` | Dimensiones, área y finitud de la superficie extraída. |

El test llamado `divergence_failure` no identifica divergencia experimental del
álabe. Comprueba un caso algebraico diseñado para ser singular.

## Tres niveles distintos

1. **Verificación de implementación:** comprobar que el código reproduce sus
   ecuaciones y referencias. Las 33 pruebas pertenecen a este nivel.
2. **Convergencia de discretización:** comprobar que una magnitud deja de
   depender materialmente del tamaño/tipo de elemento, malla VLM y paso temporal.
   No está lograda para la respuesta estructural del STEP.
3. **Validación física:** comparar predicciones con ensayos de la misma pieza,
   impresión, montaje y flujo, considerando incertidumbre. Está pendiente para
   el álabe. Se compara por separado una familia publicada de placas PLA.

La corrida histórica del STEP muestra un cambio de aproximadamente 66,1% en
desplazamiento máximo entre las dos mallas publicadas. El volumen se aproxima
mucho más al CAD, pero no permite inferir convergencia mecánica.

## Cómo ampliar el estudio de convergencia

Mantener la misma definición física de la superficie de carga, apoyos y material.
Comparar desplazamientos en puntos definidos, energía, reacciones, modos y
tensiones lejos de singularidades. Registrar la huella de la transferencia:
mantener un número fijo de vecinos puede reducir el área cargada al refinar.
Separar sensibilidad estructural de sensibilidad de la malla aerodinámica.

En transitorios, reducir el paso temporal y contrastar fase, amplitud y energía
en una ventana común. El carácter estable de Newmark para un sistema lineal
no garantiza que un paso grande resuelva las frecuencias de interés.

## Automatización

La [configuración de CI](../.github/workflows/verify.yml) ejecuta dos trabajos:

- Integridad de archivos, hashes, JSON y enlaces locales con Python estándar.
- Verificación en MATLAB R2026a sobre Ubuntu 24.04 y ejemplo FEM paramétrico, con reportes pequeños
  publicados como artefactos de Actions; no se agregan resultados al repositorio.

La CI no regenera CAD ni ejecuta siempre la malla refinada completa. La prueba
geométrica usa los CSV y la auditoría incluida. La disponibilidad de MATLAB
en Actions corresponde al mecanismo oficial para repositorios públicos de
[MathWorks](https://github.com/matlab-actions/setup-matlab#licensing).
No se necesitan secretos para este flujo público.

## Placa y respuesta de la cola: diez controles nuevos

| Control | Comparación independiente |
| --- | --- |
| `plate_rigid_modes` | Tres movimientos rígidos, sin modos espurios, masa positiva en cuadrilátero distorsionado. |
| `thin_cantilever` | Flexión de viga bajo presión uniforme, nu=0 y h/L de 1e-2 a 1e-4; verifica ausencia de bloqueo severo. |
| `simply_supported_plate` | Primera frecuencia frente a solución de Navier. |
| `pressure_patch_conservation` | Mallas no coincidentes: fuerza, primeros momentos y trabajo afín del campo de presión. |
| `harmonic_oscillator` | Amplitud y fase complejas exactas, incluido DC y resonancia amortiguada. |
| `forced_time_convergence` | Oscilador forzado exacto; convergencia temporal de segundo orden. |
| `pressure_history_bounds` | Interpolación de canales y rechazo de extrapolación. |
| `modal_solid_adapter` | Frecuencias de un sólido frente a su ruta modal original. |
| `published_data_audit` | Nueve geometrías, 27 registros, identificación del primer B y seis medias parciales documentadas. |
| `observable_comparison` | Sesgo/RMSE conocidos, normalización física, RMS y frecuencia impuesta de 5 Hz. |

`run_tail_benchmarks()` ejecuta **por separado** la comparación experimental.
Usa refinaciones [8,16] por defecto y publica errores de los tres primeros
rangos de frecuencia. Solo el primer modo B tiene identificación inequívoca
en el extracto. La columna de cambio de malla corresponde al primer modo,
no a convergencia de todas las frecuencias o tensiones.

La CI ejecuta además los dos ejemplos de cola; conserva CSV, PNG y alcance.
La comparación de placas no convierte las 33 verificaciones numéricas en
validación de la geometría STEP. Ver [auditoría y resultados](case-catalog.md).
