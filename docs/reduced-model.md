# Modelo reducido cuasiestacionario

Este módulo recibe matrices reducidas y una base cinemática. Los ejemplos
incorporados son ilustrativos; no representan parámetros identificados del CAD.
Para calcular masa y rigidez desde elementos usar el [FEM sólido](solid-fem.md).

## Ejecutar

Desde la carpeta del proyecto, en MATLAB:

```matlab
addpath('matlab')
run('matlab/inspect_reference.m') % CAD y envolventes exteriores
run('matlab/run_demo.m')          % ejemplo paramétrico completo
run('matlab/run_reference.m')     % geometría CAD + estructura ILUSTRATIVA
addpath('matlab/tests')
results = run_tests();
```

Los ejemplos guardan figuras, tablas CSV y resultados MAT en `results/`.
`run_demo` incluye un transitorio y un diagnóstico local cuasiestacionario.
**Los parámetros de ambos ejemplos son demostrativos. No son predicciones
calibradas del perfil anexado.** La procedencia estructural queda guardada
en los resultados y en los títulos de las figuras.

## Qué calcula esta versión

- Superficie aerodinámica media 3D paramétrica o interpolada desde secciones CAD.
- Aerodinámica inviscida de envergadura finita mediante VLM de herraduras.
- Masa, amortiguamiento y rigidez reducidos, incluyendo acoplamientos entre coordenadas.
- Restauración lineal más cúbica, o ley personalizada con su matriz tangente.
- Equilibrio estático con actualización geométrica y Newton con búsqueda de paso.
- Barridos de velocidad e incidencia con continuación.
- Transitorios mediante Newmark implícito y cargas **cuasiestacionarias**.
- Linealización local del mismo modelo cuasiestacionario.
- Interpolación RBF de modos muestreados y transferencia de cargas por trabajo virtual.

## Alcance físico

Esta versión **no es UVLM ni CFD**. No incluye memoria de estela, masa añadida
del fluido, separación, transición, pérdida dinámica ni resistencia viscosa.
`CDi` representa una estimación de resistencia inducida; no resistencia total.
`pressureJump` es una carga equivalente de superficie delgada, no la presión
resuelta sobre las caras del sólido.

Los transitorios sirven para variaciones suficientemente lentas como para
justificar la aerodinámica cuasiestacionaria. Evaluar la frecuencia reducida
`k = omega*c/(2*U)` forma parte de esa justificación; el código no asume un
umbral universal. **Los autovalores de `pae.linearize` no establecen por sí
solos una velocidad de flutter validada.** Para ese objetivo se requiere
un solver aerodinámico no estacionario con estados de memoria, o CFD-FSI.

La actualización de la superficie incluye cambios de normales y de influencia
aerodinámica. El modelo estructural incorporado mantiene masa constante y una
base de desplazamientos reducida; una fuerza cúbica no lo convierte en un FEM
geométricamente exacto. Las grandes rotaciones y restricciones deben modelarse
mediante `geometryFcn` y una ley estructural consistente, o un solver externo.
Una cinemática personalizada no agrega automáticamente fuerzas giroscópicas
ni masa dependiente de la configuración; esos casos requieren extender las
ecuaciones de movimiento y el integrador.
Los mecanismos con contacto, fricción seca, holgura o topes no se identifican
automáticamente. El ejemplo de deslizador es una coordenada interna acoplada
mediante rigideces: no representa la cinemática exacta del CAD.

## Formular un perfil nuevo

El modelo requiere una superficie media y datos estructurales identificados.
No asigna materiales ni apoyos a partir del nombre del archivo.

```matlab
p.spanStations = [-0.3; 0; 0.3];                 % m, orden creciente
p.leadingEdge = [0,-0.3,0; 0,0,0; 0.05,0.3,0]; % x,y,z [m]
p.chord = [0.20;0.25;0.18];                       % m
p.twist = deg2rad([0;0;-2]);                      % positivo: nariz arriba
p.camberFcn = @(xi,eta) 0.02*4*xi.*(1-xi);         % z/c
surface = pae.make_surface(p,6,20);

% Datos de FEM reducido o de identificación experimental:
data = load('mi_estructura_reducida.mat');
% data.structure contiene M,C,K,Phi,labels,source y opcionalmente cubic.
model = pae.make_model(surface,data.structure);
flow = struct('rho',1.225,'speed',12,'alpha',deg2rad(4));
equilibrium = pae.solve_static(model,flow);
[table,solutions] = pae.sweep(model,5:5:25,deg2rad([0,3,6]),flow);
```

La superficie suministrada representa **toda la envergadura finita**. Esta
implementación no aplica automáticamente simetría en la raíz, paredes de túnel
ni placas laterales. Los apoyos se incorporan a la base y matrices estructurales.
Todos los datos deben usar unidades compatibles del SI. Los ángulos de entrada
son radianes y las tablas reportan la incidencia en grados.

## Datos estructurales y grados de libertad

`structure.M`, `C` y `K` son matrices `r x r`. `M` debe ser definida positiva;
las matrices incorporadas deben ser simétricas y `C`, `K` semidefinidas positivas.
`Phi` tiene tamaño `3*N x r`, donde `N` es el número de nodos aerodinámicos.
Su orden es el de `X(:)` en MATLAB: primero todos los desplazamientos x,
después y y después z. Su escala depende de la definición de cada coordenada.
El campo opcional `coordinateUnits` registra las unidades de cada coordenada
(por ejemplo `{'m','rad'}`). Las figuras usan milímetros únicamente cuando
todas las coordenadas se han declarado en metros; en otro caso conservan la
escala original y muestran las unidades en la leyenda.

```matlab
X0 = reshape(surface.nodes,[],3);
X = reshape(X0(:) + structure.Phi*q,size(X0));
```

Una columna puede representar curvatura, torsión, flexión o el efecto de un
mecanismo. Una coordenada interna del deslizador puede tener columna nula en
`Phi` y aun influir mediante entradas no diagonales de `K`, `M` o una ley propia.
Para bloquear un grado de libertad se deben reducir consistentemente matrices
y base; no usar una rigidez numéricamente infinita.

Si los modos de un FEM se conocen en otra malla:

```matlab
targetPoints = reshape(surface.nodes,[],3);
Phi = pae.interpolate_basis(samplePoints,sampleModes,targetPoints);
% samplePoints: Nsample x 3; sampleModes: Nsample x 3 x r.
% Conservar las mismas coordenadas, normalización y matrices modales M,C,K.
```

La interpolación usa RBF cúbica con reproducción afín y hasta 1500 muestras
representativas. No elimina discontinuidades físicas ni identifica uniones;
muestrear regiones compatibles y revisar la interpolación. Para datos planos,
los puntos objetivo deben pertenecer al mismo subespacio afín.

## Ecuaciones y acoplamiento

El sistema reducido resuelve:

```text
M qddot + C qdot + fint(q) = Qa(q,qdot;U,alpha) + Qexterna
fint(q) = K q + cubic .* q.^3                     [ley incorporada]
X = X(q),   J = dX(:)/dq
Qa = J' * f_nodal,   f_nodal = Hforce' * F_panel
```

El operador adjunto conserva trabajo virtual para el mapa elegido. Las fuerzas
de Kutta-Joukowski se aplican en el punto medio del vórtice ligado a un cuarto de
la cuerda de cada panel. La condición de no penetración se impone a tres cuartos.
La estela de cada herradura es recta, larga y alineada con la corriente libre.
No evoluciona con el tiempo. `coreRadius` y `wakeLength` son opciones del VLM.

`solve_static` verifica el residual de equilibrio; si no converge, falla de
forma explícita. Puede devolver el diagnóstico sin lanzar error con
`struct('errorOnFailure',false)`. En un barrido, las filas no convergidas quedan
marcadas. No interpretar sus coeficientes como equilibrios válidos.

## Extender los solvers

```matlab
structure.geometryFcn = @(q) mi_geometria(q); % devuelve [X,J]
structure.internalFcn = @(q) mi_estructura(q); % devuelve [fint,Kt]
structure.aeroFcn = @(nodes,V,flow) mi_aerodinamica(nodes,V,flow);
```

El callback aerodinámico devuelve fuerzas por panel y campos `CL,CDi,Cm,fidelity`
para el posprocesamiento. La interfaz actual es de cargas instantáneas. Un
UVLM/CFD transitorio **requiere extender la integración para almacenar y
actualizar estados del flujo**: no basta sustituir el callback por una función
con estado oculto que cambia en cada iteración de Newton.

## STEP suministrado y preparación de otro CAD

El CAD original no se distribuye en el repositorio; los ejemplos consumen los
CSV preparados. La lectura
directa mediante `fegeometry` de MATLAB R2026a falló con
`Unable to triangulate one or more faces`, también al reducir la desviación
relativa. Se preparó una triangulación mediante OpenCASCADE/Gmsh. El solver
MATLAB consume el caché ya incluido y no necesita Python para ejecutarlo.

El STEP contiene un sólido BREP con 587 caras. **Eso no especifica qué regiones
deben estar pegadas, deslizar, contactar o tener materiales diferentes.** La
triangulación sirve para inspección geométrica, no sustituye una malla estructural
con propiedades y restricciones.

Para preparar un STEP nuevo, en un entorno Python con `numpy` y `gmsh`:

```bash
python -m pip install numpy gmsh
python tools/prepare_step.py mi_perfil.step data/mi_perfil --axes 0,2,1 --mesh-size 0.006
```

`--axes` indica los índices de los ejes CAD de cuerda, envergadura y altura
respectivamente; 0=X, 1=Y, 2=Z. Gmsh convierte la geometría a metros. Las secciones
se extraen intersectando triángulos con planos, y luego tomando las envolventes
superior e inferior. No se emplean los puntos de control del STEP como contorno.
Este procedimiento requiere secciones representables por curvas superiores e
inferiores monovaluadas; revisar accesorios, cuerpos externos y discontinuidades.
No se identifican automáticamente borde de ataque, sentido de flujo ni ejes.
En el CAD de referencia se observan pequeñas discontinuidades de la envolvente
inferior cerca de la zona interna del deslizador, especialmente en las secciones
cercanas a los extremos. Revisar la figura `cad_inspection.png` y definir la piel
aerodinámica nominal antes de usar una malla fina como representación física.

```matlab
[surface,info] = pae.surface_from_sections('data/mi_perfil',6,20,"constant");
```

La política `constant` extiende las secciones de 1% y 99% a los extremos de la
envergadura. Es una aproximación explícita que hay que revisar para alas no
prismáticas o puntas especiales. `sampled` analiza solamente el tramo muestreado.
El cambio de ejes y origen usado por la superficie debe aplicarse también a los
modos estructurales. La cuerda incluye toda la envolvente exterior: un accesorio
que sobresalga podría invalidar la identificación aerodinámica.

## Comprobaciones y validación posterior

La suite verifica simetría, límite de sustentación de gran alargamiento, escala
con densidad/velocidad, refinamiento de malla, trabajo virtual con una geometría
rotacional, equilibrio y linealización contra un caso analítico, rechazo de un
problema estático singular, oscilador amortiguado exacto, conservación de un
equilibrio en el transitorio, reproducción afín de modos y escala del CAD.

Estos son controles del código y sus aproximaciones. Para validar el perfil
real se necesitan ensayos carga-deformación, modos/amortiguamiento, condiciones
de apoyo y comparación de fuerzas y deformaciones en túnel. También se requiere
contrastar el VLM con aerodinámica viscosa en el régimen de Reynolds pertinente.

## Referencias metodológicas

- [NASA: Optimum lattice arrangement developed from a rigorous analytical basis](https://ntrs.nasa.gov/citations/19760021094).
- [SHARPy: arquitectura de acoplamiento aeroelástico no lineal](https://doi.org/10.21105/joss.01885). El paquete propio `uvlm` incorpora esa familia de ecuaciones con
  alcance y verificación propios; no usa SHARPy como motor ni reproduce todas sus capacidades.
- [Ochi et al. (2024): validación de un perfil pasivo](https://doi.org/10.2514/1.J063000).
- [MathWorks: geometría FEM](https://www.mathworks.com/help/pde/ug/fegeometry.html) y [análisis estructural lineal](https://www.mathworks.com/help/pde/structural-mechanics.html).
