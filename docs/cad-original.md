# Campaña estructural del CAD original

El objeto principal es **Airfoil with slider alone.stp**, conservando su
conectividad. El [informe actual](cad-original-results.md) documenta el cierre
numérico bajo hipótesis autorizadas de material, apoyo y guía. Las placas/vigas
son controles independientes. El [primer bloque](cad-original-history.md) y sus
mallas fallidas siguen disponibles como historia, sin sustituir los resultados actuales.

## Reproducir las tres mallas principales

Instalar las dependencias opcionales de `requirements-preprocessing.txt` en un
entorno Python. Gmsh/OpenCASCADE solo leen y mallan el CAD; el FEM es propio.
Comprobar el SHA del STEP contra el contrato antes de regenerar:

```sh
python3 tools/mesh_step_extruded.py "/ruta/Airfoil with slider alone.stp" \
  results/cad_original/meshes/extruded_2 --size .004 --min-size .0001 --layers 14 \
  --volume-audit data/reference/cad_volume_audit.json
python3 tools/mesh_step_extruded.py "/ruta/Airfoil with slider alone.stp" \
  results/cad_original/meshes/extruded_3 --size .0025 --min-size .000075 --layers 20 \
  --volume-audit data/reference/cad_volume_audit.json
python3 tools/mesh_step_extruded.py "/ruta/Airfoil with slider alone.stp" \
  results/cad_original/meshes/extruded_4 --size .0018 --min-size .00006 --layers 22 \
  --volume-audit data/reference/cad_volume_audit.json
```

El preprocesador comprueba caps planos y extrusión, conserva todas las curvas
originales y divide los prismas de manera conforme. No cura ni elimina detalles
del BREP. Las razones jacobianas muy pequeñas del CAD se conservan y se informa
su efecto en los [residuos numéricos](numerical-residuals.md).

```matlab
addpath('matlab')
run_verification();
study=run_cad_original('cases/cad_original/extruded-fine.json',false);
run_cad_mechanisms(study.output);
run_cad_guided_modes(study.output,1,3);
run_cad_guided_modes(study.output,2,3);
run_cad_guided_modes(study.output,3,3);
run_cad_modal_shape_comparison(study.output);
run_cad_full_reduced(fullfile(study.output,'tet10_level_3'));
```

Cada campaña usa un directorio nuevo con UTC, contrato, hashes y copia exacta
de fuentes. El nivel mayor tiene 1 134 945 DOFs y requiere factorización dispersa;
se ejecutó en un equipo de 24 GiB, en serie para las grandes matrices. El
presupuesto distingue matrices, tripletas por bloque y llenado de Cholesky.
Las mallas de resultados se regeneran desde el STEP; no se publican silenciosamente
como si estuvieran incluidas entre las dos mallas Tet4 históricas.

La configuración principal es [extruded-fine.json](../cases/cad_original/extruded-fine.json).
`run_cad_original` prepara A; `run_cad_mechanisms` añade las hipótesis B/C por
restricciones medias ponderadas por área sobre la barra 7/9. No separa nodos
ni inventa una articulación de la conectividad. El
[contrato de hipótesis](../cases/cad_original/assumptions.json) especifica su alcance.

## Carga, regiones y sensores

Los ejes del solver son X,Z,Y del CAD. La permutación refleja la orientación:
las caras se orientan de nuevo desde el tetraedro propietario y los momentos
se transforman con `det(R)*R`, respecto al mismo origen. El caso usa material
isótropo; una base ortótropa exige transformar la constitutiva y verificarla.

Las caras CAD 568–585 y 403 son candidatas a piel superior posterior; 407–425,
a piel inferior. Estas etiquetas se vinculan al hash del STEP y a firmas de
área, centroide y caja envolvente. No asignan automáticamente función mecánica
a nervaduras, flexuras o slider. El parche de carga se limita a
`0.78 <= x/c <= 0.84`, `0.25 <= eta <= 0.75`.

El operador recorta los triángulos por ese rectángulo físico, integra en cada
subtriángulo y evalúa las funciones de forma de la cara original. Verifica
cobertura proyectada, normal superior y exposición mediante rayos verticales
en una retícula fija; las paredes de cavidades no se cargan por cercanía.
Las presiones positivas actúan hacia dentro. El área superficial real y las
fuerzas horizontales pueden variar con la aproximación geométrica; el área
proyectada y la fuerza vertical son idénticas entre mallas. No se renormaliza
una fuerza para ocultar un error de área.

`H` lleva desplazamientos nodales a puntos de integración. `H'` lleva las
fuerzas de esos puntos a nodos y conserva trabajo virtual; se comprueban
fuerza y momento con las coordenadas físicas. La comparación con 16 vecinos
utiliza las mismas fuerzas puntuales y solo nodos de la piel seleccionada.
El informe contiene radio de huella y cambio de los nueve sensores.

Los sensores están en `x/c=0.999`, nueve estaciones `eta=0.05..0.95`, en
material expuesto de la piel superior. Ese desplazamiento está **0,1% de cuerda
antes del borde geométrico**, para evitar la ambigüedad de una arista singular.
Se informa `w_z` firmado. La posición x/y y el operador de observación se fijan
antes de deformar, con interpolación de la cara; no se siguen nodos distintos
por cercanía. Su altura inicial varía con la aproximación geométrica de cada
malla. Los sensores Tri6 curvos usan inversión paramétrica Newton con
comprobación de pertenencia y residuo físico.

## Tet10 y límites geométricos

`fem3d.elevate` añade un nodo compartido por arista: vértices 1–4, seguidos
por aristas 12,23,31,14,24,34. No cambia el STEP ni proyecta automáticamente
nodos al CAD. El piloto compara Tet4/Tet10 sobre la misma geometría facetada.

`tet10` implementa gradientes, deformaciones, rigidez, masa consistente y
volumen. En geometría afín usa rigidez exacta de grado 2 y masa exacta de grado
4. En geometría curva emplea integración Duffy de pesos positivos, con orden
refinable, orientación consistente y controles adicionales del jacobiano.
La verificación incluye orden de integración 4/6/8; cada geometría curva nueva
requiere su propio estudio. Los controles discretos no prueban positividad
continua del jacobiano en cualquier geometría arbitraria.

Ensamblaje por bloques, caras Tri6, restricciones de todos los nodos,
conectividad, tensiones de centroide y visualización están adaptados. La
visualización subdivide Tri6 para dibujar; esa subdivisión no altera la carga
ni la rigidez. Las tensiones son de centroide del elemento y no máximos nodales
ni resistencia validada del filamento.

## Comparación geométrica curva

```sh
python3 tools/project_step_quadratic.py "/ruta/Airfoil with slider alone.stp" \
  results/cad_original/meshes/extruded_4 results/cad_original/meshes/curved_projection
```

```matlab
run_cad_geometry_comparison(fullfile(study.output,'tet10_level_3'), ...
    'results/cad_original/meshes/curved_projection');
```

La proyección coloca nodos intermedios en curvas/caras originales y rechaza
ambigüedades o inversión en los controles. La comparación utiliza una nueva
solución FEM completa y la misma presión de 100 Pa. No se usa una base afín
congelada como medida del error geométrico.

`run_cad_curved_quadrature(curvedFolder)` compara órdenes Duffy 4 y 6 en
todos los elementos, usando los seis modos y el desplazamiento estático de
esa solución curva. La tolerancia es 0,1% para energías, masas y operadores
proyectados escalados. Es una comprobación de integración de campos físicos;
no sustituye una nueva solución completa a orden 6 ni acota todos los vectores
del espacio FEM. Se conserva además el control geométrico de jacobianos.

## Resultados y aceptación

Cada nivel guarda atlas, sensores, curva presión–deformación 0/25/50/100 Pa,
reacciones, energías, modos, FRF, matrices/estados y presupuesto de recursos.
`source_hashes.csv`, `input_hashes.csv` y `executed_sources` trazan lo ejecutado.
Las tensiones de centroide no son una resistencia validada del filamento.

Se aceptan por separado el cambio entre mallas, el error algebraico en energía,
la carga transferida y la banda modal. La sexta frecuencia del modelo A no
cumple 1% y queda fuera de la banda principal 0–60 Hz. `structuralModelClosed`
se refiere a cierre numérico condicional; la identificación experimental es
otra cuestión. Los indicadores históricos originales se preservan y el
[informe](cad-original-results.md) explica su interpretación actual.

La [ruta acoplada](cad-coupled.md) parte de K/M y sensores de este mismo CAD.
La [extensión viscosa](viscous-flow.md) declara por separado Reynolds,
representación de cavidades y error de discretización. Para revisar archivos
aún sin seguimiento Git: `python3 tools/check_repository.py --include-working-tree`.


`run_cad_posterior(campaign,dynamicFolder,basisFile)` extrae líneas materiales
posteriores, desplazamientos, pendiente y curvatura ajustada en ventanas físicas.
Conserva por separado la convergencia de esa derivada y la de TE/energía; la
curvatura ajustada no se interpreta como tensión en una arista singular.
