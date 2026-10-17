# Solver FEM propio y referencia histórica desacoplada

Esta versión calcula la estructura a partir de **nodos, elementos y propiedades
del material**. `+fem3d` implementa el elemento sólido, el ensamblaje, los apoyos,
el equilibrio, la recuperación de tensiones, los modos y el transitorio.
`+pae` aporta el VLM de herraduras ya implementado en el proyecto.
No se necesita introducir matrices reducidas elegidas a mano.

El núcleo numérico usa MATLAB base: álgebra lineal dispersa y operadores
geométricos. No llama al solver FEM de PDE Toolbox, ANSYS ni a un solver CFD.
Para el STEP, Gmsh/OpenCASCADE se utiliza **exclusivamente para leer el CAD y
generar conectividades**; el mallado paramétrico de los ejemplos también tiene
una implementación propia en MATLAB. Los CSV permiten ejecutar una malla
preparada sin tener Gmsh o Python instalados.

**Estado del caso STEP:** las dos mallas Tet4 descritas en esta guía son la
referencia histórica que no convergió. La [campaña actual del CAD](cad-original-results.md)
usa tres niveles Tet10 y cumple los criterios declarados de cola, energía y
banda fundamental. Incluye guías, reducción, acoplamiento estático/dinámico y
extensión viscosa. Las propiedades y el montaje siguen siendo hipótesis.
La verificación del código y la precisión de una discretización son cuestiones distintas.

## Ejecutar y examinar el código

Desde la carpeta del proyecto:

```matlab
addpath('matlab')
run('matlab/run_fem_demo.m')       % sólido NACA paramétrico, apoyos supuestos
run('matlab/run_fem_reference.m')  % STEP suministrado, sólido y VLM de una vía
run('matlab/run_fem_convergence.m')% dos mallas STEP; puede tardar
addpath('matlab/tests')
tests = run_fem_tests();          % verificaciones del FEM y de la transferencia
```

El caso del archivo anexado se ejecuta con `run_fem_reference.m`. Revisar
`data/fem_reference/mesh_info.json`: la malla tiene volumen auditado y elementos
muy delgados que requieren convergencia estructural. El ejemplo es un
**sólido homogéneo unido**, con el 30% anterior de la
cuerda empotrado. Ambos son supuestos de desarrollo. El archivo STEP no define
los grados de libertad, contacto, fricción ni apoyos del deslizador.

Los resultados se guardan en `results/fem_demo/` y `results/fem_reference/`.
Las figuras indican la amplificación del dibujo; los colores representan
desplazamientos reales del modelo, sin esa amplificación. Los archivos MAT
contienen malla, matrices y solución para poder inspeccionar cada operación.

La auditoría independiente `data/reference/cad_volume_audit.json` confirma un
BREP válido y un volumen de `6.155285293e-4 m³`. La malla inicial difiere un
0,476% y la refinada un 0,102% respecto a la integración adaptativa del CAD.
El cálculo no adaptativo de OpenCASCADE/Gmsh daba un volumen aproximadamente
12% menor que el adaptativo; no era una referencia adecuada para este sólido.
Se conserva esa diferencia y el método utilizado en los metadatos. La
[documentación de BRepGProp](https://dev.opencascade.org/doc/refman/html/class_b_rep_g_prop.html)
describe ambas rutas de integración. La auditoría de volumen comprueba la
geometría; no sustituye la convergencia de desplazamientos y tensiones.

Con el mismo material, incidencia de 3°, velocidad de 10 m/s y apoyo supuesto:

| Malla | Nodos / tetraedros | Desplazamiento máximo [mm] | Energía elástica [J] |
|---|---|---:|---:|
| Inicial | 21.149 / 69.082 | 0,005520 | 3,124e-7 |
| Refinada | 85.367 / 286.037 | 0,016286 | 9,767e-7 |

La diferencia de desplazamiento relativa a la malla refinada es **66,1%**.
Las masas y volúmenes casi coinciden, pero la respuesta mecánica no. Se necesita
mejorar los elementos en las zonas delgadas y revisar la transferencia local,
además de identificar el apoyo y el mecanismo reales. El máximo de tensión
también cambia mucho; no se utiliza para certificar resistencia.

Para generar otra malla STEP, usar un entorno Python con `gmsh` y `numpy`.
La auditoría independiente necesita opcionalmente `cadquery-ocp` en un
intérprete compatible (la auditoría suministrada se ejecutó con Python 3.12):

```bash
python tools/audit_step_volume.py mi_alabe.step data/mi_perfil/cad_volume_audit.json
python tools/mesh_step_volume.py mi_alabe.step data/mi_malla \
  --size 0.006 --min-size 0.0002 --volume-audit data/mi_perfil/cad_volume_audit.json
```

El CAD original no se publica. El código no modifica el STEP propio que se
proporcione como entrada. Los datos del mallado preparado ya
están incluidos; esas dependencias solo se necesitan al regenerarlo.

## Qué significa desacoplado aquí

La secuencia es:

1. Calcular las fuerzas aerodinámicas sobre la geometría inicial a una velocidad
   e incidencia especificadas.
2. Transferir esas fuerzas a la malla estructural conservando fuerza, momento
   y trabajo virtual.
3. Resolver el equilibrio estructural con esas cargas fijas.

La aerodinámica no se recalcula sobre la forma deformada en este caso. Es una
respuesta estructural a cargas aerodinámicas, o FSI de una vía. Constituye una
base verificable para el acoplamiento aeroelástico posterior; **no calcula
divergencia ni flutter**.

Los modos y el Newmark de `+fem3d` son estructurales. El transitorio acepta una
historia de fuerzas externa, sin memoria del fluido. El VLM es inviscido y de
superficie delgada: entrega resultantes y un salto de presión equivalente,
no los dos campos de presión de las caras reales. Para el régimen de Reynolds
de un túnel pequeño, la transición y la separación requieren validación
aerodinámica adicional.

## Formulación del FEM

Cada tetraedro tiene cuatro nodos, tres desplazamientos por nodo y doce grados
de libertad. El orden es:

```text
ue = [ux1 uy1 uz1 ux2 uy2 uz2 ux3 uy3 uz3 ux4 uy4 uz4]'
epsilon = [exx eyy ezz gamma_xy gamma_yz gamma_xz]'
sigma   = [sxx syy szz tau_xy tau_yz tau_xz]'
```

Las funciones de forma son `N=[1-r-s-t, r, s, t]`. La matriz jacobiana tiene
como columnas las aristas `X2-X1`, `X3-X1` y `X4-X1`. Se obtienen las derivadas
espaciales de las funciones de forma y se construye `B`, de modo que
`epsilon = B*ue`. El elemento lineal tiene deformación constante.

```text
Ve = abs(det(J))/6
Ke = Ve * B' * D * B
Me = rho*Ve/20 * kron(ones(4)+eye(4), eye(3))
```

La masa es consistente, con integración analítica. La rigidez también se
integra exactamente para este elemento y un material constante. Las
contribuciones se ensamblan en matrices dispersas globales, por bloques para
limitar memoria. `tet4.m` contiene el elemento y `assemble.m` el ensamblaje.

`isotropic.m` obtiene `D` de `E` y `nu`. `orthotropic.m` invierte una matriz de
flexibilidad simétrica construida con nueve constantes independientes. Sus
ejes materiales coinciden con los ejes globales; una orientación de impresión
distinta exige transformar el tensor o suministrar `D` ya transformada. Se
comprueba que `D` sea definida positiva. Cada elemento puede tener su propio
`materialId`, lo que permite asignar piel, nervaduras y zonas flexibles.

Los desplazamientos prescritos se eliminan exactamente:

```text
Kff * uf = Ff - Kfc * uc
Rc = (K*u - F)c
```

No se usan rigideces artificialmente grandes para empotrar nodos. Se verifica
que cada sólido conectado tenga restringidos sus seis movimientos rígidos y
que la matriz libre admita una factorización de Cholesky. Se reportan residual,
balance de fuerzas, balance de momentos y energía elástica.
La factorización dispersa aplica una permutación AMD para reducir el relleno;
se recupera la solución en el orden original de los grados de libertad. Esto
es relevante para las mallas grandes del STEP.

`stress.m` calcula `epsilon_e=B*ue` y `sigma_e=D*epsilon_e`. Se conserva la
tensión por elemento, sin un suavizado que oculte discontinuidades. Von Mises
es una medida para visualizar tensiones; **no es un criterio de fallo validado
para piezas FDM**. No se calcula un factor de seguridad con resistencia última
como si fuera un límite elástico.

El análisis modal resuelve `Kff*phi = omega^2*Mff*phi` y normaliza cada modo a
masa unitaria. Se reporta el residual del problema propio.
`newmark.m` integra `M*uddot+C*udot+K*u=F(t)` con `beta=1/4`, `gamma=1/2`.
Por defecto `C=0`: la literatura consultada no aporta amortiguamiento para tu
pieza. Se pueden introducir dos coeficientes de Rayleigh identificados.

## Transferencia de cargas

El VLM calcula fuerzas en puntos del vórtice de cuarto de cuerda. `transfer.m`
construye una interpolación local `H` que reproduce desplazamientos afines:

```text
u_aero = H * U_structural
F_structural = H' * F_aero
H*ones = ones
H*X_structural = X_force_points
```

De estas identidades se obtiene conservación de fuerza y momento; el uso de la
traspuesta conserva trabajo virtual. El procedimiento se comprueba
independientemente en los tests. Puede producir pesos negativos, por lo que
debe inspeccionarse el tamaño y la selección de cada parche.

Esta transferencia conserva resultantes pero no identifica la piel física ni
reconstruye presión CFD. Para un mecanismo con piezas móviles distintas, se
deben limitar los nodos candidatos a la pieza que recibe cada carga. La
selección automática de nodos próximos podría atravesar una holgura o cargar
una nervadura interna incorrectamente. `traction.m` permite aplicar
directamente presión o tracción sobre caras del sólido cuando se dispone de
ese campo: presión positiva actúa hacia dentro sobre normales exteriores.
Las fuerzas concentradas equivalentes del VLM tampoco permiten interpretar
sin más los picos de tensión locales: al refinar cambia su distribución nodal.
Para estudiar concentraciones de tensión físicas, usar cargas superficiales
distribuidas y una transferencia cuya huella represente el área cargada.

## Materiales de manufactura aditiva

Las referencias incorporadas proceden de la tabla 3 de
[Malim et al. (2023)](https://doi.org/10.3390/aerospace10030256), un estudio de
palas impresas que identifica propiedades en probetas y compara deformaciones
de piezas representativas. La configuración corresponde a FDM con relleno
100%, raster −45/+45° y capa de 0,2 mm.

| Referencia del código | E [GPa] | rho [kg/m³] | Resistencia última X/Y/Z [MPa] |
|---|---:|---:|---|
| `tough_pla_malim2023` | 2,617 | 1160 | 47 / 47 / 38 |
| `abs_malim2023` | 1,813 | 1084 | 37,5 / 37,5 / 20 |

En ambos casos se adopta `nu=0.35` como **supuesto explícito**: no se presenta
como medición de ese estudio. Las resistencias se guardan como metadatos.
La aproximación isotrópica es una referencia inicial para esa configuración,
no una ley universal para cualquier pieza impresa.

Si se cambia la orientación, el relleno, el filamento o la temperatura, se
requiere identificar propiedades efectivas. No se multiplica automáticamente
`E` por el porcentaje de relleno. Para TPU o resinas elastoméricas tampoco se
infiere `E` de Shore A: hacen falta curvas tensión-deformación y, para grandes
deformaciones, un modelo constitutivo apropiado.

## Literatura que orienta el desarrollo

Esta es una revisión enfocada de fuentes primarias, no una revisión sistemática
exhaustiva. Los trabajos de morphing pasivo orientan la cinemática y validación;
los de impresión aportan propiedades y sus límites de transferencia.

| Trabajo | Evidencia pertinente | Consecuencia para este solver |
|---|---|---|
| [Taguchi et al., 2020](https://doi.org/10.2322/tjsass.63.18) | Perfil pasivo con zona flexible desde `x/c=0.3`; medida óptica de deformación, fuerzas y presión en túnel. | Contrastar simultáneamente la forma y las cargas. El apoyo anterior del ejemplo es una idealización inspirada en esta distribución, no la réplica de su montaje. |
| [Kai et al., 2023; disponible online desde 2022](https://doi.org/10.1007/978-981-19-2689-1_78) | NACA0012 pasivo con mecanismo compliant impreso y evaluación aerodinámica en túnel. | La flexibilidad del mecanismo interno debe modelarse junto con la piel; sustituirla por un ala maciza cambia la respuesta. |
| [Ochi et al., 2024](https://doi.org/10.2514/1.J063000) | Simulación y validación de un perfil pasivo impreso; el resumen describe acoplamiento FEM–fluido y comparación Euler/RANS. | La aerodinámica viscosa es una referencia para validar las cargas del VLM. No se extrajeron parámetros estructurales de este trabajo. |
| [Dogea, Yan y Millar, 2023](https://doi.org/10.1007/s44245-023-00013-8) | Nervadura NACA0018 impresa en PLA, ensayos y FEM. La respuesta depende de apoyos, posición de carga y recortes. | Los apoyos y la arquitectura interna requieren tanta atención como el valor de E. |
| [Gonabadi, Yadav y Bull, 2020](https://doi.org/10.1007/s00170-020-06138-4) | Ensayos de PLA FFF muestran dependencia de orientación y relleno, con comportamiento anisótropo. | Incorporar ortotropía cuando las propiedades direccionales se conozcan; no trasladar indiscriminadamente propiedades de filamento. |
| [Malim et al., 2023](https://doi.org/10.3390/aerospace10030256) | Palas impresas y cálculo FSI estacionario de una vía. | Referencia de material para iniciar el caso desacoplado. |
| [Huang y Zhang, 2026](https://doi.org/10.1016/j.ast.2026.112449) | Comparación modal experimental y numérica de modelos viga, placa y mixtos para morphing de camber. | Elegir la idealización estructural según flexión en cuerda, envergadura y torsión; una viga única no reproduce todos esos efectos. |

Además, [el conjunto de datos de Huang y Zhang, versión 2, junio de
2026](https://data.mendeley.com/datasets/8jxw5533py/2) publica ensayos de PLA
impreso en siete orientaciones, datos modales, tolerancias, CAD y archivos FEM.
La ruta nueva incorpora un extracto de nueve geometrías de placa y compara
frecuencias con una formulación MITC4 propia. Ver [casos publicados](case-catalog.md).
Esta comparación no resuelve la convergencia del sólido STEP ni identifica su slider.

La conclusión de implementación es que MATLAB resulta adecuado para construir
y verificar el solver propio. La dificultad principal es elegir una cinemática
y discretización que representen el mecanismo, resolver paredes muy delgadas y
validar propiedades impresas; no el ensamblaje disperso por sí solo.

## Definir un caso general

```matlab
mesh.nodes = X;                 % N x 3, metros
mesh.tetra = T;                 % Ne x 4, índices empezando en 1
mesh.materialId = ids;          % Ne x 1
[mesh.faces,mesh.faceOwners] = fem3d.boundary(T,X);
material = fem3d.material_reference("tough_pla_malim2023");
system = fem3d.assemble(mesh,material);

% Empotramiento de un grupo de nodos:
bc.dofs = fem3d.dofs(nodos_fijos,1:3);
bc.values = zeros(size(bc.dofs));

% Presión constante sobre las caras elegidas:
F = fem3d.traction(mesh,caras_cargadas,100);   % Pa
solution = fem3d.solve_static(system,F,bc);

% O aerodinámica desacoplada, con una superficie definida mediante +pae:
flow = struct('rho',1.225,'speed',10,'alpha',deg2rad(3));
study = fem3d.one_way(system,surface,flow,bc,nodos_de_la_piel);
```

Para una guía ideal alineada con x, restringir y/z en los nodos apropiados con
`fem3d.dofs(nodos_guia,[2,3])`. Esto solo prescribe desplazamientos: no crea
automáticamente contacto entre cuerpos, fricción ni una articulación del CAD.
Para sólidos independientes deben existir nodos y conectividades independientes
antes de añadir las restricciones físicas entre ellos.

## Verificación y límites de precisión

`run_fem_tests.m` comprueba movimientos rígidos, masa y energía cinética,
invariancia de orientación, patch test afín con tracciones analíticas, equilibrio
de presión, refinamiento de una viga en voladizo, modos, Newmark, transferencia
conservativa y el límite isotrópico de la ley ortótropa. También rechaza apoyos
insuficientes y audita la malla de referencia. Los CSV conservan los errores,
además del resultado de cada test.

En la verificación de flexión, el error frente a Euler–Bernoulli baja al refinar
la malla. Esto **no garantiza precisión en cualquier piel delgada**: los Tet4
pueden ser excesivamente rígidos en flexión y en materiales casi incompresibles.
Una solución con residual pequeño satisface el sistema discreto, pero todavía
puede tener un gran error de discretización. Deben compararse mallas usando
desplazamiento, energía, modos y tensiones en zonas relevantes, no solamente
el máximo puntual de tensión junto a un empotramiento.

## Extensiones necesarias para un estudio aeroelástico completo

1. **Discretización para flexibilidad:** tetraedros cuadráticos para sólidos o
   elementos de lámina/viga para pieles y ligamentos. Para láminas de Mindlin,
   una formulación como MITC evita el bloqueo por cortante; el
   [curso de Bathe sobre placas y láminas](https://ocw.mit.edu/courses/res-2-002-finite-element-procedures-for-solids-and-structures-spring-2010/resources/lecture-20/)
   explica esta dificultad. El espesor de pared debe representarse explícitamente.
2. **Mecánica del mecanismo:** contactos, guías y restricciones multipunto a
   partir de la arquitectura real. El STEP unido no basta para identificarlas.
3. **Grandes deformaciones:** equilibrio incremental con una cinemática
   consistente, fuerzas internas y matriz tangente. Usar un material muy blando
   en este FEM lineal no añade esa capacidad.
4. **Acoplamiento estático de dos vías:** actualizar la superficie, recalcular
   cargas y converger simultáneamente desplazamiento y residual de fuerzas,
   con relajación cuando corresponda.
5. **Aeroelasticidad dinámica:** aerodinámica no estacionaria con memoria
   (UVLM o CFD), amortiguamiento identificado y estudio de estabilidad.

La primera etapa implementada permite separar y verificar la estructura, la
aerodinámica y la transferencia antes de introducir estas extensiones.

## Extensión al CAD original y Tet10

La [campaña CAD](cad-original.md) añade Tet10, caras Tri6, integración física
de presión y sensores materiales de la cola. Los ejemplos anteriores mantienen
su configuración; la campaña nueva guarda resultados en otra carpeta.
La verificación incorpora flexión y modos analíticos de viga, masa,
orientación, geometría curva y conservación de transferencia. La elevación
de aristas no mejora por sí sola los elementos esbeltos ni demuestra
convergencia del álabe. Consultar el estado de ejecución del
[plan CAD](plan-cad-original.md) antes de usar sus valores como predicción.
