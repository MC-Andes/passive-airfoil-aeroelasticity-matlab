# Datos, geometría y preprocesamiento

## Qué se incluye y por qué

Los archivos de `data/` son entradas de ejecución. Permiten correr todos los
ejemplos MATLAB sin un kernel CAD ni el archivo original de Inventor. Los
CSV son texto, índices basados en uno y coordenadas en metros.

| Ruta | Consumidor |
| --- | --- |
| `reference/cad_info.json` y `section_01.csv` a `section_05.csv` | Superficie VLM y tests de geometría. |
| `reference/cad_nodes.csv`, `cad_triangles.csv` | Inspección geométrica. |
| `reference/cad_volume_audit.json` | Referencia independiente de volumen y procedencia. |
| `fem_reference/{nodes,tetra,faces}.csv`, `mesh_info.json` | Caso sólido inicial y test de auditoría. |
| `fem_refined/{nodes,tetra,faces}.csv`, `mesh_info.json` | Diagnóstico de convergencia. |
| `manifest.json` | Hashes y dimensiones de las entradas publicadas. |

No se publican el STEP/Inventor original, STL redundantes, segmentos auxiliares
de intersección ni resultados. El hash del STEP utilizado se conserva como
procedencia: `5bacb470df1b1239a0d8af2181f10ca14077a4c0eba3b3531126c8981ea4793e`.
Ese hash no permite recuperar el CAD ni demuestra por sí solo acceso al archivo.
La auditoría incluida documenta una ejecución previa sobre esa fuente.

## Formato de malla volumétrica

| Archivo | Encabezado y significado |
| --- | --- |
| `nodes.csv` | `x_m,y_m,z_m`: coordenadas en ejes CAD. El número de fila define el ID del nodo. |
| `tetra.csv` | `n1,n2,n3,n4,cad_volume`: conectividad y etiqueta de volumen CAD. |
| `faces.csv` | `n1,n2,n3,cad_face`: triángulos de superficie y etiqueta CAD. |
| `mesh_info.json` | Unidades, cantidades, método, tamaño objetivo y auditoría de volumen. |

`read_mesh` conserva las etiquetas CAD, vuelve a extraer caras desde los
tetraedros y asigna `materialId=1` inicialmente. Un volumen CAD no constituye
automáticamente una región material. Debe existir compatibilidad de índices,+nodos y conectividades; no se admiten nodos huérfanos en el ensamblaje.

`fem3d.read_mesh(folder,[1,3,2],surface.canonicalOriginCAD)` aplica la convención
del ejemplo. Para un CAD diferente hay que identificar sus ejes y origen.
El lector usa el volumen relativo registrado y rechaza una diferencia de 1%
o más. Ese control no sustituye un análisis de calidad y convergencia.

## Secciones para la aerodinámica

`section_XX.csv` contiene `xi,x_m,upper_m,lower_m,camber_m`. `xi` es la
coordenada normalizada de cuerda; las otras columnas están en m en el plano CAD
de cuerda/altura. `camber_m=(upper_m+lower_m)/2`.

Los cortes se obtuvieron al 1%, 25,5%, 50%, 74,5% y 99% de la envergadura.
Las envolventes descartan huecos interiores y suponen caras superior/inferior
monovaluadas. Un accesorio o una discontinuidad puede invalidar esa lectura
aerodinámica aunque el sólido CAD sea válido.

`pae.surface_from_sections` obliga a elegir `constant` o `sampled`. La primera
extiende las secciones próximas a los extremos hasta la envergadura física;
la segunda usa solo el intervalo muestreado. No hay una identificación
automática fiable de puntas o del contorno nominal.

## Auditoría independiente de volumen

El BREP de referencia fue aceptado por `BRepCheck_Analyzer`. La integración
adaptativa de OpenCASCADE produjo `6.155285293484151e-4 m³`. El cálculo no
adaptativo produjo `5.430025132465182e-4 m³`; la diferencia numérica del kernel
era importante para estas superficies spline.

| Malla | Volumen [m³] | Diferencia relativa respecto del adaptativo |
| --- | ---: | ---: |
| Inicial | 6,125988658e-4 | 0,00475959 |
| Refinada | 6,148991689e-4 | 0,00102247 |

El volumen geométrico mejora al refinar, pero el desplazamiento aún no
converge. Ambas comprobaciones deben presentarse por separado.

## Preparar otro CAD

Crear un entorno Python 3.12 separado e instalar los requisitos opcionales:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-preprocessing.txt
```

En Windows, activar con `.venv\Scripts\activate` desde la terminal apropiada.
El sistema debe soportar las bibliotecas nativas de Gmsh/OpenCASCADE; las
herramientas son opcionales y no se usan en la CI del solver.

Preparar envolventes, auditar y generar volumen para un archivo proporcionado
por el usuario:

```bash
python tools/prepare_step.py mi_alabe.step data/mi_perfil --axes 0,2,1 --mesh-size 0.006
python tools/audit_step_volume.py mi_alabe.step data/mi_perfil/cad_volume_audit.json
python tools/mesh_step_volume.py mi_alabe.step data/mi_solido \
  --size 0.006 --min-size 0.0002 \
  --volume-audit data/mi_perfil/cad_volume_audit.json
```

El preprocesador convierte el CAD a metros. El auditor de OpenCASCADE fija su
unidad interna a milímetro y convierte el volumen a m³. El auditor y el mallador
deben operar sobre exactamente el mismo archivo para que coincida el SHA-256.

`prepare_step.py` también puede generar STL y segmentos auxiliares para
diagnóstico; esos archivos no son necesarios para ejecutar el solver y no se
versionan en la referencia publicada. Las dependencias nativas conservan sus
propias licencias y no se redistribuyen dentro del repositorio.

## Reproducibilidad de las dos mallas históricas

Las mallas incluidas son entradas congeladas y verificadas por hash. Sus
metadatos describen Gmsh 4.15.2. La malla inicial se produjo con tamaño máximo
0,015 m, mínimo 0,0015 m, curvatura 8 y optimización Netgen. La refinada usó
0,006 m, 0,0002 m, curvatura 16 y sin ese paso adicional. La herramienta actual
utiliza curvatura 16, por lo que sus opciones de tamaño por sí solas no
reproducen exactamente la triangulación histórica inicial. Las versiones del
kernel y el algoritmo también pueden cambiar numeración y conectividades.

Para una reproducción idéntica utilizar los CSV publicados. Para un nuevo
mallado guardar las opciones completas y repetir auditoría y convergencia,
en lugar de esperar igualdad binaria con las mallas históricas.

## Integridad de los datos

`tools/check_repository.py` compara cada entrada con `data/manifest.json` y
verifica tamaños/tablas. Al modificar una malla, actualizar el manifiesto,
documentar el método y ejecutar MATLAB otra vez. Un hash de integridad detecta
cambios de archivo; no valida la física, la calidad del mallado ni la asignación
de materiales y condiciones de frontera.
