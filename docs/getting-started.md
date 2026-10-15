# Instalación y ejecución

## Requisitos

- MATLAB R2026a es la versión verificada localmente. La CI solicita R2026a.
- El cálculo usa MATLAB base. No requiere PDE Toolbox, Optimization Toolbox,
  Simulink ni un solver FEM externo.
- Las mallas y secciones numéricas necesarias están en `data/`.
- Python y las bibliotecas CAD son opcionales: solo se usan para regenerar
  entradas a partir de un archivo STEP propio.

No se garantiza compatibilidad con versiones anteriores a R2026a sin ejecutarlas.
Los scripts usan, entre otras funciones, `arguments`, `readmatrix`,
`exportgraphics`, tablas y álgebra dispersa. No es un proyecto verificado en Octave.

## Obtener el proyecto

```bash
git clone https://github.com/MC-Andes/passive-airfoil-aeroelasticity-matlab.git
cd passive-airfoil-aeroelasticity-matlab
```

No hacen falta Git LFS, submódulos ni descargas de datos adicionales. Los datos
CSV son entradas de referencia; los resultados se generan localmente.

## Primera ejecución en MATLAB

Cambiar la carpeta actual a la raíz del clon. Después:

```matlab
addpath('matlab')
report = run_verification();
run('matlab/run_fem_demo.m')
```

La primera llamada ejecuta las 23 verificaciones y crea reportes bajo
`results/verification/`. La segunda calcula un caso paramétrico, sus primeros
modos y una vibración estructural libre. Las figuras del FEM se exportan sin
abrir una ventana interactiva.

En una terminal donde `matlab` esté en `PATH`:

```bash
matlab -batch "addpath('matlab'); run_verification(); run('matlab/run_fem_demo.m');"
```

Si MATLAB no está en `PATH`, usar la ruta al ejecutable de su instalación o
abrir MATLAB desde su interfaz. No es necesario agregar rutas personales al
código del proyecto.

## Caso basado en el STEP de referencia

```matlab
run('matlab/inspect_reference.m')
run('matlab/run_fem_reference.m')
run('matlab/run_fem_convergence.m')
```

La inspección consume la triangulación incluida. El FEM consume la malla
volumétrica preparada, no necesita el STEP original. El tercer script compara
dos mallas y muestra que el caso todavía no converge. Puede requerir más tiempo
y memoria que el ejemplo paramétrico. No hay una estimación universal de RAM:
la factorización dispersa depende de la conectividad y los apoyos.

Los ejemplos son scripts: crean variables en el workspace y guardan las
salidas con nombres fijos. Al repetirlos se reemplazan sus reportes anteriores.
Guardar aparte una corrida que se quiera conservar.

## Comprobar la publicación sin MATLAB

```bash
python3 tools/check_repository.py
```

Este control usa solo la biblioteca estándar de Python. Comprueba el contenido
versionado, hashes de entradas, metadatos, enlaces locales de documentación y
sintaxis de los scripts Python. No sustituye los tests de MATLAB ni regenera
las mallas.

## Qué archivos se publican

Se incluyen código MATLAB, pruebas, herramientas opcionales de preparación,
documentación, metadatos y CSV necesarios para los ejemplos. Se excluyen
`results/`, `.cache/`, entornos, archivos MAT/FIG, ZIP, miniaturas, el CAD
original y exports redundantes. No se incluye ninguna copia de MATLAB, Gmsh,
OpenCASCADE ni de artículos científicos.

## Problemas habituales

| Síntoma | Comprobación |
| --- | --- |
| `Unable to resolve the name fem3d...` o `pae...` | Ejecutar `addpath('matlab')` desde la raíz; no agregar `+fem3d` o `+pae` directamente. |
| Falta `nodes.csv` o `cad_info.json` | Verificar que se clonó todo `data/`, conservando los subdirectorios. |
| `fem3d:Unconstrained` | Revisar apoyos, cuerpos desconectados y movimientos rígidos restantes. |
| `fem3d:MeshVolume` | Revisar unidades, auditoría y tolerancia de volumen; no forzar la aceptación. |
| El campo de tensiones depende mucho de la malla | Evaluar discretización, picos en cargas/apoyos y distribución física de presión. |
| `ModuleNotFoundError` en una herramienta CAD | Instalar las dependencias en el entorno Python opcional, no en MATLAB. |
| Falla la lectura STEP | Los ejemplos preparados pueden ejecutarse; para otro CAD revisar cierre, unidades y tolerancias en el preprocesador. |

Para reportar un fallo, adjuntar comando, versión, commit y mensaje exacto en
un issue. Los reportes deben distinguir un error de ejecución de una limitación
física del modelo.
