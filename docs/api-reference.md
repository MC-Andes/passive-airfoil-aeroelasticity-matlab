# Referencia de API

Todas las funciones son MATLAB base. Los índices de nodos y componentes son
enteros positivos basados en uno. Las magnitudes emplean SI. Los argumentos
entre corchetes en las descripciones son opcionales; no representan sintaxis
literal de MATLAB. `help fem3d.tet4`, por ejemplo, muestra el contrato junto al
código implementado.

## Estructuras principales

| Estructura | Campos y dimensiones |
| --- | --- |
| `mesh` | `nodes`: N×3; `tetra`: Ne×4; `materialId`: Ne×1 opcional; `faces`: Nf×3 y `faceOwners`: Nf×1 para cargar/representar fronteras. |
| `material` | `D`: 6×6 simétrica definida positiva; `rho`: escalar positivo. Se pueden agregar `name`, `E`, `nu`, procedencia y parámetros de impresión. |
| `system` | `K`,`M`: 3N×3N dispersas; `mesh`, `materials`, `materialId`, `volume`: Ne×1, `totalMass`, `dofOrder`, `formulation`. |
| `bc` | `dofs`: vector único de grados restringidos; `values`: vector del mismo tamaño, cero si se omite. |
| `surface` | `nodes`: (nx+1)×(ny+1)×3; `referenceArea`, `referenceChord`, `span`; mallas paramétricas incluyen `xi`,`eta`. |
| `flow` | `rho` del aire, `speed` no negativa, `alpha` rad; opcionales `beta`, `referencePoint`, `referenceArea`, `referenceChord`. |

La densidad del material y la densidad del aire son entradas diferentes.

## Geometría y materiales del sólido

| Función | Uso y retorno |
| --- | --- |
| `fem3d.box_mesh(lengths,divisions)` | Caja de longitudes `[Lx,Ly,Lz]`, con subdivisiones `[nx,ny,nz]`. Origen en cero; devuelve `mesh`. |
| `fem3d.airfoil_mesh(chord,span,thicknessRatio,divisions)` | NACA00xx macizo, envergadura centrada; divisiones por defecto `[30,16,4]`. Trunca x/c a 0,005 y 0,995 para evitar celdas degeneradas. |
| `fem3d.structured_mesh(grid)` | `grid` de ni×nj×nk×3, con al menos dos nodos por dirección. Divide cada hexaedro lógico en seis tetraedros conformes. |
| `fem3d.read_mesh(folder,axisOrder,origin)` | Lee cuatro archivos del formato publicado. Por defecto `[1,2,3]` y `[0,0,0]`. Recalcula la orientación de las caras de frontera. |
| `[faces,owners]=fem3d.boundary(T,X)` | Extrae caras exteriores orientadas y su elemento propietario. Rechaza caras compartidas por más de dos tetraedros. |
| `fem3d.isotropic(E,nu)` | Devuelve `D`; requiere E>0 y −1<nu<0,5. |
| `fem3d.orthotropic(E,nu,G)` | `E=[Ex,Ey,Ez]`, `nu=[nu_xy,nu_yz,nu_xz]`, `G=[Gxy,Gyz,Gxz]`. Los ejes materiales coinciden con los globales. |
| `fem3d.material_reference(name)` | `tough_pla_malim2023` por defecto o `abs_malim2023`; propiedades y procedencia, con nu supuesto. |

`nu_ij` significa deformación transversal en j bajo tracción uniaxial en i.
La simetría de la flexibilidad impone los coeficientes recíprocos. No se incluye
una rotación automática de la ortotropía.

## Elemento, ensamblaje y restricciones

| Función | Uso y retorno |
| --- | --- |
| `[Ke,Me,B,V]=fem3d.tet4(X,D,rho)` | X de 4×3, `Ke`,`Me` de 12×12, `B` de 6×12 y volumen positivo. Tolera orientación invertida; rechaza degeneración. |
| `fem3d.assemble(mesh,materials)` | Ensambla K/M por bloques de 3000 elementos. Un arreglo de materiales permite seleccionar por `materialId`. Rechaza nodos huérfanos y propiedades inválidas. |
| `fem3d.dofs(nodes,components)` | Devuelve índices intercalados; componentes 1=x, 2=y, 3=z, por defecto 1:3. |
| `fem3d.prepare_constraints(nd,bc)` | Valida restricciones y agrega `bc.free`. Las restricciones son desplazamientos constantes. |
| `fem3d.check_supports(mesh,bc)` | Comprueba seis movimientos rígidos restringidos en cada componente conectado. No detecta todos los mecanismos internos. |

Para usar dos materiales, crear estructuras con campos compatibles:

```matlab
m1 = struct('D',fem3d.isotropic(2.617e9,.35),'rho',1160);
m2 = struct('D',fem3d.isotropic(1.813e9,.35),'rho',1084);
mesh.materialId = ones(size(mesh.tetra,1),1);
mesh.materialId(elementos_region_2) = 2;
system = fem3d.assemble(mesh,[m1,m2]);
```

`elementos_region_2` se obtiene de una partición física definida por el usuario.
No equivale automáticamente a una lista de volúmenes CAD.

## Cargas y solución estructural

### `F=fem3d.traction(mesh,faces,load)`

Devuelve fuerzas nodales N×3. `load` puede ser presión escalar constante,
presiones Nf×1, tracción global 1×3 o tracciones Nf×3, siempre en Pa. Para una
sola cara, un vector 1×3 significa tracción. Las normales deben estar orientadas
hacia fuera para que presión positiva corresponda a carga hacia dentro.

### `solution=fem3d.solve_static(system,F,bc)`

Acepta F de N×3 o un vector 3N×1 intercalado. Elimina restricciones y usa
Cholesky disperso con ordenamiento AMD. Devuelve:

- `u`: vector intercalado y `displacement`: N×3.
- `F`: carga aplicada en el orden intercalado; `reaction`: N×3, cero en grados libres.
- `bc`, `freeResidual`, `relativeResidual`, `forceBalance`, `momentBalance`.
- `strainEnergy`, `maximumDisplacement`, `factorNonzeros` y `stress`.

`relativeResidual` se normaliza con el máximo del lado derecho libre y `eps`.
No usarlo como estimador del error de discretización.

### `stresses=fem3d.stress(system,u)`

Devuelve `strain`, `sigma` de Ne×6, `vonMises` de Ne×1,
`maximumVonMises` y el orden de componentes. No suaviza entre elementos.

### `modes=fem3d.modal(system,bc,count)`

Seis modos por defecto. Requiere apoyos homogéneos y menos modos que grados
libres. Devuelve `frequencyHz`, `omega`, `Phi` de 3N×count,
`eigenResidual` y `bc`. Cada modo cumple masa unitaria. La fase/signo del modo
se fija con su componente de mayor magnitud; dos modos cercanos pueden cambiar
de orden o base entre plataformas sin modificar su subespacio físico.

### `response=fem3d.newmark(system,bc,time,forceFcn,u0,v0,rayleigh)`

`time` debe ser creciente. El callback `forceFcn(t)` devuelve N×3 o 3N×1.
`u0`, `v0` son vectores intercalados; por defecto, reposo compatible con los
apoyos. `rayleigh=[aM,bK]` es no negativo, por defecto `[0,0]`.
Se usa `C=aM*M+bK*K`, con aM en 1/s y bK en s. No hay amortiguamiento identificado
para el álabe real.

Devuelve `time`, `U`,`V`,`A` de 3N×Nt, energía estructural, residuales libres y
parámetros de integración. Aquí `U` es una historia vectorial, no la matriz
N×3 de la solución estática. Condiciones prescritas variables con el tiempo,
masa dependiente de configuración y no linealidad no están implementadas.

## Transferencia y ejecución desacoplada

### `[H,info]=fem3d.transfer(nodes,points,candidates,neighbors)`

Nodos estructurales N×3 y puntos de carga P×3 en el mismo sistema de ejes y
origen. Todos los nodos son candidatos por defecto; `neighbors=16`. Amplía
parches si no reproduce un movimiento afín, hasta los límites del código.
H es P×N dispersa. `H*U` interpola desplazamientos y `H'*forces` transfiere
fuerzas. `info` registra defecto afín, radios de parche y pesos negativos.

Seleccionar candidatos de la piel/pieza correcta. La conservación de
resultantes no impide un parche físicamente incorrecto ni elimina picos
originados por cargas concentradas.

### `study=fem3d.one_way(system,surface,flow,bc,candidates)`

Los candidatos por defecto son nodos de `system.mesh.faces`. Usa los valores
de área y cuerda de `surface` como referencia de coeficientes. Devuelve
`aero`, `structure`, `H`, `transferInfo`, defectos de fuerza/momento/trabajo,
`flow` y descripción de acoplamiento. No calcula modos ni transitorios: se
solicitan por separado.

### `fig=fem3d.plot_result(system,solution,titleText,scale)`

Representa desplazamiento y tensiones exteriores, con amplificación geométrica
`scale=1` por defecto. Los valores de color no se amplifican. Crea una figura
invisible que se puede exportar con `exportgraphics`.

## Aerodinámica común: paquete `pae`

| Función | Contrato principal |
| --- | --- |
| `make_surface(planform,nx,ny)` | Planform con `spanStations`, `leadingEdge` Nstation×3 y `chord`; opcionales `twist`, `camberFcn`. Por defecto nx=6, ny=16. |
| `surface_from_sections(folder,nx,ny,edgePolicy)` | Lee envolventes; `edgePolicy` obligatorio: `constant` o `sampled`. Devuelve superficie e info del CAD. |
| `panels(nodes)` | Normales, área, vórtices ligados A/B, control y fuerzas, conectividades y operadores `Hforce`,`Hcontrol`. |
| `segment_velocity(points,a,b,core)` | Velocidad inducida de un segmento recto de circulación unitaria, regularizado. |
| `vlm(nodes,nodeVelocity,flow,options)` | Velocidades nodales N×3. Opciones de vórtice: `coreRadius`, `wakeLength`. Devuelve fuerzas, circulación, CL, CDi, Cm y residuales. |

La presión equivalente es `pressureJump`, no Cp de las caras del sólido. A
velocidad cero las fuerzas son cero y los coeficientes normalizados son NaN.
El momento se toma respecto de `flow.referencePoint` o `[0,0,0]` si se omite;
`Cm` usa su componente y. El VLM no devuelve resistencia viscosa.

## Modelo reducido: paquete `pae`

| Función | Contrato principal |
| --- | --- |
| `make_model(surface,structure)` | Requiere `M,C,K,Phi,labels,source`; masa definida positiva y estructura compatible. |
| `geometry(model,q)` | Devuelve nodos deformados y Jacobiano de la cinemática. |
| `internal_force(model,q)` | Fuerza interna y tangente; ley lineal/cúbica o callback. |
| `map_loads(J,Hforce,panelForces)` | Fuerzas generalizadas por trabajo virtual y cargas nodales. |
| `evaluate(model,q,qdot,flow)` | Estado coherente de geometría, cargas y fuerzas internas para una evaluación. |
| `solve_static(model,flow,q0,overrides)` | Equilibrio con Newton, tangente aerodinámica y búsqueda de paso; falla si no converge por defecto. |
| `sweep(model,speeds,angles,flowBase,overrides)` | Continuación; devuelve tabla y soluciones. Las incidencias de entrada están en radianes. |
| `linearize(model,solution,overrides)` | Matriz local de estados y autovalores cuasiestacionarios; no es flutter validado. |
| `simulate(model,flowInput,time,q0,v0,overrides)` | Newmark con masa constante y cargas cuasiestacionarias. `flowInput` es estructura o función del tiempo. |
| `jacobian(fun,x,step)` | Diferencias centrales; `step` obligatorio, escalar o uno por coordenada. Usa step·max(1,abs(x)) por componente. |
| `interpolate_basis(samplePoints,sampleModes,targetPoints)` | RBF cúbica afín, 3 a 1500 muestras; devuelve Phi en el orden `X(:)`. |
| `options(overrides)` | Valida opciones; valores por defecto en la tabla siguiente. |
| `example_structure(surface)` | Estructura ilustrativa de cuatro coordenadas, no identificada del CAD. |
| `demo_case(nx,ny)` | Geometría y modelo ilustrativo; por defecto nx=4, ny=12. |
| `plot_study(model,summary,solutions,response)` | Gráficos del barrido y respuesta opcional; esta interfaz pertenece al modelo reducido. |

| Opción reducida | Valor |
| --- | --- |
| `maxIterations` | 35 |
| `relativeTolerance` | 1e-8 |
| `absoluteTolerance` | 1e-9 |
| `differenceStep` | 1e-6 |
| `maxLineSearch` | 12 |
| `errorOnFailure` | true |
| `newmarkBeta`, `newmarkGamma` | 0,25 y 0,5 |

Estas opciones no son parámetros del FEM lineal `fem3d.solve_static`.

## Entrada de verificación

`report=run_verification()` ejecuta ambas suites y escribe CSV/JSON bajo
`results/verification/`. Un test fallido interrumpe la ejecución con error.
La interfaz no descarga dependencias ni necesita una sesión gráfica.
