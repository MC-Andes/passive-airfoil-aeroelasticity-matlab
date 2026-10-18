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

`report=run_verification()` ejecuta las cuatro suites y escribe CSV/JSON bajo
`results/verification/`. Un test fallido interrumpe la ejecución con error.
La interfaz no descarga dependencias ni necesita una sesión gráfica.

## Placa y adaptador desacoplado de presión

La formulación, convenciones y ejemplos se detallan en
[borde de salida](trailing-edge.md) y [casos publicados](case-catalog.md).
Estas funciones usan MATLAB base y apoyos homogéneos.

| Función | Contrato |
| --- | --- |
| `plate.mesh(x,y,t)` | Vectores crecientes en metros; espesor escalar o callback de centroides; devuelve cuadriláteros antihorarios. |
| `plate.shape(xi,eta)` | N y derivadas paramétricas bilineales; nodos (-1,-1),(1,-1),(1,1),(-1,1). |
| `plate.mitc4(X,t,E,nu,rho)` | K12×12, M12×12, cargas por Pa12×1 y área; X4×2. DOFs intercalados w,theta_x,theta_y. |
| `plate.assemble(mesh,material)` | K, M dispersas, P por presión elemental, masa y áreas. Material E,nu,rho. |
| `plate.pressure_map(mesh,xEdges,yEdges)` | B3N×Nparches para presiones rectangulares; x más rápido. Fuerza y primeros momentos conservados para ese campo. |
| `plate.probe(mesh,points,component)` | O para puntos Nsensores×2; componente 1=w, 2=theta_x, 3=theta_y. |
| `ael.static(K,fixed,F,O)` | Equilibrio completo; F puede tener varias columnas. Devuelve U,Y,reacciones y residual relativo. |
| `ael.reduce(K,M,fixed,count,zeta,B,O)` | Modos normalizados en masa, matrices modales M/C/K/B/O, Phi, frecuencias y residuos. Sin suponer orden físico de DOFs. |
| `ael.harmonic(model,fHz,input)` | Amplitudes pico complejas; input Ncanales×1 o ×Nfrecuencias. Exp(i*omega*t). |
| `ael.transient(model,time,inputFcn,q0,v0)` | Newmark beta=1/4,gamma=1/2; estados iniciales modales opcionales cero. Q/V/A,Y,energía,residuo. |
| `ael.pressure_history(time,pressure)` | Callback desde presiones Ncanales×Ntiempos Pa; interpolación lineal, sin extrapolación. |
| `ael.response_metrics(time,Y)` | Señales Nsensores×Ntiempos uniformes: media, RMS fluctuante, pico a pico, espectro Hann y bin dominante. |
| `ael.compare_series(t,Y,tRef,Yref,scale)` | Interpola simulación a observaciones cubiertas; bias,RMSE,RMSE/scale,error máximo. Sin ajuste de fase. |
| `ael.reference_cases()` | Nueve structs con geometría/material nominal, datos experimentales y procedencia. |
| `ael.build_plate_case(c,nChord)` | Sistema placa y DOFs de apoyo para el contrato de caso publicado; incorpora interfaz de espesores. |

`ael` no sustituye las funciones cuasiestacionarias de `pae`; integra una
estructura lineal bajo entradas prescritas. B/O de sólidos deben construirse
con su orden intercalado de traslaciones; los de placa contienen pendientes.

## Sólido cuadrático y campaña CAD

La ruta completa y sus criterios se describen en [campaña CAD](cad-original.md).

| Función | Contrato |
| --- | --- |
| `fem3d.elevate(mesh)` | Eleva Tet4 a Tet10 conforme, con aristas rectas; no proyecta al CAD. Conserva etiquetas de caras. |
| `fem3d.tet10(X,D,rho,order)` | K30×30, M30×30, B de centroide, volumen e informe de jacobianos. Orden Duffy por eje por defecto 4; geometría afín sin argumento de orden usa integrales exactas. |
| `fem3d.tet10_shape(r)` | Funciones cuadráticas y derivadas; aristas 12,23,31,14,24,34. |
| `fem3d.tri_shape(r,count)` | Tri3/Tri6; aristas de cara 12,23,31. |
| `fem3d.surface_quadrature(mesh,faces,order)` | H, puntos, normales y pesos de área en caras explícitas Tri3/Tri6, incluida geometría curva. |
| `fem3d.traction(mesh,faces,load)` | Presión positiva hacia dentro o tracción cartesiana. También acepta callback `@(points,normals)` que devuelve Nq×1 o Nq×3. |
| `fem3d.plot_faces(faces)` | Triángulos de visualización y mapa a caras físicas; Tri6 se divide en cuatro. |
| `cad.load_case(root,path)` | Contrato versionado; valida auditoría de fuente y escenario A; calcula Re y b/c. |
| `cad.transform(axisOrder,origin)` | Matriz R de vectores polares y det(R)R de momentos axiales; origen expresado en ejes permutados. |
| `cad.face_labels(mesh)` | Une caras FEM a etiquetas CAD por conectividad, rechazando pérdidas o duplicados. |
| `cad.pressure_patch(mesh,config)` | Recorte físico de piel exterior Tri3/Tri6; integración, H y fuerza por Pa, exposición y balances; delega geometría curva en `cad.curved_patch`. |
| `cad.sensors(mesh,config)` | H y O de nueve puntos materiales cerca del TE, sin ajuste a nodo cercano; inversión paramétrica Newton sobre Tri6 curvos. |
| `cad.atlas(mesh,config,output)` | Firmas por cara CAD, localización de elementos esbeltos e intervalos verticales de material. |
| `cad.resource_estimate(mesh)` | Aristas reales, DOFs Tet10 y cota de tripletas sin ensamblar el sistema. |
| `cad.modal_compare(a,b,eta)` | MAC ponderado y asignación de modos observados; marca cercanía modal que exige subespacios. |
| `run_cad_original(casePath,runTet10)` | Campaña sin sobrescrituras: tres niveles del elemento configurado, piloto opcional, curva de presión, modos, FRF condicional y aceptación. |

`assemble`, `boundary`, `read_mesh`, `check_supports`, `stress` y `plot_result`
aceptan ahora Tet4/Tet10 y Tri3/Tri6. El CSV de tetraedros contiene 4 o 10
nodos y una etiqueta de volumen; el de caras, 3 o 6 nodos y una etiqueta CAD.
La recuperación de tensiones de Tet10 es de centroide. `solve_static` usa
Cholesky disperso equilibrado y refinamiento con el operador original; guarda
su historia de residuos. `modal` equilibra K/M mediante la misma congruencia
y calcula residuos en las coordenadas originales.

`run_cad_modal_response(solutionFolder)` calcula un decaimiento modal libre
a partir de un `solution.mat` CAD, compara dos pasos temporales con la
solución analítica amortiguada y guarda señales, RMS/espectro, energía y
alcance. El amortiguamiento es supuesto y la perturbación inicial impuesta;
no simula una excitación del flujo ni valida la reducción espacial.


## Mecánica, FSI y estela

| Interfaz | Contrato |
| --- | --- |
| `fem3d.solve_constrained`, `solve_contact`, `solve_friction` | Restricciones exactas, conjunto activo unilateral y fricción incremental con disipación. |
| `fem3d.nonlinear_element`, `solve_nonlinear` | Sólido total Lagrangiano Saint Venant–Kirchhoff, tangente consistente y continuación. |
| `fem3d.modal_constrained`, `newmark_observed` | Modos guiados sin penalización y observación de transitorio completo sin guardar todos los estados. |
| `cad.build_reduced`, `remap_reduced` | Base del FEM CAD con flexibilidad residual; nueva malla fluida conserva K/M y observación. |
| `fsi.static`, `dynamic`, `linear_dynamic` | Equilibrio con geometría actualizada, partición fuerte UVLM/Newmark y tangente con memoria completa. |
| `uvlm.initialize`, `step`, `linear_model` | Anillos/Kutta, estela congelada o libre, cargas estacionarias/no estacionarias y tangente. |
| `uvlm.wake_vertex_gradient` | Derivada analítica del vértice compartido por dos segmentos; contraste complejo independiente. |
| `fsi.compact_dynamic`, `truncate` | Eliminación Schur exacta y subconjunto modal consistente de todos los operadores. |
| `nsflow.project_penalized`, `external` | Proyección variable presión/Brinkman y flujo exterior conservativo 2D/3D. |
| `nsflow.surface_traction` | Presión y esfuerzo viscoso 2D fuera de la máscara, con referencia Cp aguas arriba explícita. |
| `run_cad_aeroelastic` | Campaña estática/dinámica con dt, paneles, núcleo, longitud de estela y alcance declarados. |
| `run_cad_viscous`, `run_cad_viscous_postprocess` | Campañas Re/malla/dominio/reinicio y perfiles físicos comunes de estela. |
| `build_native_kernels` | Compilación opcional de kernels C propios equivalentes; todo el solver admite alternativa MATLAB. |

Las [guías acoplada](cad-coupled.md) y [viscosa](viscous-flow.md) describen
unidades, hipótesis, conservación y reproducción. Los runners de sensibilidad
no son calibraciones del espécimen ni modelos de turbulencia.


## Estudio actualizado PETG y servo

| Función | Contrato |
| --- | --- |
| `morphing.load_case(root,[path])` | Schema 2, STEP actualizado y datos de estudio coherentes; no reutiliza el loader histórico. |
| `morphing.boundary_conditions(mesh,c,command_m)` | Apoyo sobre cuatro caras de alojamiento hasta inserción declarada; desplazamiento uniforme de interfaz del slider. Rechaza SHA/firmas incompatibles. |
| `morphing.observers(mesh,c)` | Puntos materiales en envergadura positiva; offset TE referido a extensión real del CAD, no a una cuerda nominal inexistente. |
| `morphing.aerodynamic_coefficients(points,forces,c,U,alpha_deg)` | Puntos [m] y fuerzas sobre cuerpo [N], canónicas; CL/CD de viento y Cm sobre el eje, S=.09 m², c=.25 m. |
| `run_morphing_study([casePath])` | Genera 231 filas de planificación; no calcula CFD ni inventa recorrido final. |
| `run_morphing_mechanics([casePath],[order])` | Prueba tangente completa, modos bloqueados y tensiones. Tet10 por defecto; Tet4 es comparación diagnóstica. |
| `run_morphing_nonlinear([command_m],[casePath])` | SVK Tet4 con continuación y Cauchy/Green; amplitud pequeña por defecto, sin calificar resistencia. |
| `run_morphing_review(manifestPath)` | Comparación de campañas explícitas, hashes y figuras; no hereda cierre del CAD anterior. |

`fem3d.nonlinear_tet4_batch` calcula la misma fuerza, tangente y energía del
elemento total Lagrangiano existente, con páginas y memoria acotadas. El flag
`mesh.nonlinearBatch=true` lo activa para Tet4. Su contraste independiente y
la solución analítica de una barra con desplazamiento finito impuesto forman
parte de los seis controles morphing. El solver normaliza también con las
reacciones de los DOFs impuestos cuando la fuerza exterior es cero.
