# Mecánica y flujo acoplado del CAD original

Esta ruta mantiene el sólido conectado del STEP y el FEM/fluido MATLAB propios.
El usuario autorizó sustituir las entradas experimentales desconocidas por
hipótesis geométricas y bibliográficas. El contrato de
[hipótesis](../cases/cad_original/assumptions.json) distingue esas entradas de
las dimensiones auditadas. No existe validación experimental de esta pieza.

## Estructura y mecanismo

`fem3d.tet10` integra el sólido isoparamétrico, con masa consistente y frontera
Tri6. `tet10_batch` implementa la misma integración por páginas para geometría
curva; se contrasta con el ensamblaje elemental independiente. La elevación
`fem3d.elevate` produce aristas rectas. La herramienta
`tools/project_step_quadratic.py` proyecta nodos de frontera sobre caras/curvas
originales, comprueba entidades compartidas y jacobianos, y conserva un informe
separado. Una proyección rechazada no sustituye silenciosamente una malla.

`cad.ray` y `cad.face_point` invierten el mapeo Tri6. El parche plano en
coordenadas físicas conserva sus dimensiones al refinar; para piel curva se
integra en x/y con cuadratura compuesta y el jacobiano de la superficie. Su
refinamiento se controla con `loads.integrationDivisions`. El test de
paraboloide contrasta fuerza y momento analíticos. En todos los casos la
interpolación H y la fuerza H' son adjuntas.

Las caras CAD 7/9 delimitan la barra interna candidata a slider. El operador
`cad.slider_constraints` ajusta sus seis coordenadas rígidas mediante
integración superficial ponderada por área, sin depender del número de nodos.
A conserva el sólido; B bloquea las seis coordenadas; C permite su traslación
media en cuerda. No se separan nodos ni se introducen bisagras. La interpretación
como guía es una hipótesis, no una identificación del ensamblaje real.

`solve_constrained` resuelve restricciones afines mediante multiplicadores.
`constraint_response` reutiliza una factorización del sólido monolítico para
su mapa de flexibilidad de guía. `solve_contact` verifica gap y reacción por
conjunto activo; `solve_friction` conserva la posición tangencial anterior y
aplica Coulomb incremental, con disipación comprobada. Son contactos de pequeñas
excursiones; no representan deslizamiento finito entre superficies desconocidas.

La formulación no lineal total Lagrangiana usa Green–Lagrange y St
Venant–Kirchhoff con tangente material/geometría consistente. Se verifican
objetividad a rotación finita, derivada direccional y barra a 5% de extensión.
Es una ley elástica supuesta; no describe plastificación o viscoelasticidad
identificadas de PLA. Los indicadores u/L y rotación local determinan cuándo
compararla con la aproximación lineal en cada campaña.

## Aerodinámica y transferencia

`cad.surface` obtiene las dos pieles expuestas por intersección, y su superficie
media material. Excluye el primer/último 0,1% de cuerda para evitar la inversión
singular del borde. El operador medio divide las fuerzas entre ambas pieles
como idealización delgada. Conserva trabajo, fuerza y momento; no resuelve Cp
superior/inferior, presión interna o cavidades.

`fsi.vlm` convierte fuerzas de herraduras aplicadas en cuarto de cuerda a
fuerzas nodales con el operador adjunto correspondiente. La presión equivalente
del panel no se interpreta como una presión uniforme con otro brazo de momento.
Las mallas aerodinámica y estructural se refinan por separado.

`fsi.static` actualiza posiciones y normales, itera cargas/estructura con
relajación y comprueba residuos de carga y desplazamiento. El fallo de iteración
se informa como tal. Se contrastan alta rigidez, carga pequeña y dos factores
de relajación antes de interpretar realimentación aeroelástica.

## UVLM y dinámica

El paquete `uvlm` calcula anillos ligados, no penetración, emisión de la
circulación TE del paso anterior y convección de estela congelada o libre.
Conserva todos los anillos emitidos y rechaza exceder el presupuesto de memoria;
no descarta circulación. Las fuerzas incluyen Joukowski de segmentos ligados
y Bernoulli no estacionario rho*A*n*dGamma/dt. El modelo no representa separación
viscosa ni una estela turbulenta.

`run_uvlm_frequency_benchmarks` obtiene la solución periódica del sistema UVLM
discreto, con sus anillos de estela, y la compara con Theodorsen y Sears. La
prueba transitoria de heave contrasta por otra ruta esa respuesta. Se usa
k=omega*c/(2U), exp(i*omega*t), y Sears referido al borde de ataque mediante
exp(-i*k). Espacio, dt, horizonte y núcleo se registran por separado.

`uvlm.linear_model` deriva la tangente del UVLM propio alrededor del equilibrio
CAD. Conserva historias de circulación y de posiciones materiales TE. El
contraste con UVLM no lineal se realiza primero en un ala de prueba y después
en movimiento prescrito de la superficie real. Esta tangente es adecuada para
perturbaciones pequeñas; no se extrapola a ciclos límite de amplitud finita.

`fsi.linear_dynamic` resuelve simultáneamente el Newmark estructural y las
cargas UVLM de cada paso. La convolución se agrupa en productos matriciales
exactos para reutilizar memoria y la factorización aerodinámica. El contraste
bloques 1/16 verifica que no hay compresión del núcleo ni pérdida de historia.
`fsi.dynamic` conserva asimismo una interfaz particionada con subiteraciones,
siempre copiando el mismo estado fluido anterior hasta converger.

La respuesta libre comienza con una perturbación TE de amplitud declarada.
No se añade presión sinusoidal para producir vibración. El ajuste de crecimiento
usa la primera coordenada modal, excluye dos picos iniciales y amplitudes
inferiores a 0,1% del máximo; guarda ventana, picos y residuo del ajuste.
La mezcla de modos en una señal TE tardía puede producir una frecuencia falsa,
por lo que esa señal se conserva pero no se utiliza ciegamente para identificar
flutter. Son necesarios suficientes ciclos y convergencia antes de clasificarlo.

## Ejecución

Desde la raíz:

```matlab
addpath('matlab')
run_verification();
% Los directorios de campaña se pasan explícitamente; nunca se elige un
% resultado por fecha sin comprobar su contrato y aceptación.
run_cad_mechanisms(campaignFolder);
run_cad_full_reduced(fullfile(campaignFolder,'tet10_level_3'));
cad.build_reduced(fullfile(campaignFolder,'tet10_level_3'),16,10,outputFolder);
run_cad_aeroelastic(fullfile(outputFolder,'reduced_cad.mat'));
run_uvlm_frequency_benchmarks();
```

La [ruta viscosa](viscous-flow.md) estudia otras magnitudes físicas. Los resultados,
convergencias y limitaciones actuales se concentran en el
[informe del CAD](cad-original-results.md). Los ejemplos de placas de literatura
continúan siendo verificaciones/comparaciones separadas del mecanismo original.


## Coste y reproducción de las refinaciones grandes

`build_native_kernels()` compila opcionalmente dos kernels C propios de
Biot–Savart y su derivada, y añade `.cache/native` a la sesión. La alternativa
MATLAB sigue disponible al retirar ese directorio del path. Los kernels
comparten ecuaciones, regularización y cortes; se contrastan por separado
contra MATLAB y derivación compleja. No delegan el fluido en una biblioteca
externa. Las campañas registran hash del binario y fuente C ejecutados.

`staticRelaxation=1` acelera los casos de realimentación débil; los residuos
siguen siendo los mismos y la dependencia de relajación se comprueba aparte.
`checkStaticLimits=false` permite repetir exclusivamente una refinación sin
recalcular los límites ya verificados; su archivo declara que esos controles
no se ejecutaron en esa repetición. La aceptación cita las campañas separadas.
No se sustituye un fallo de convergencia por un resultado parcialmente iterado.


## Revisión conjunta de campañas

`cad.compare_free_response(coarseFolder,fineFolder,speed)` compara las mismas
señales físicas con interpolación temporal: cambio de frecuencia/decaimiento,
error máximo normalizado por el pico de referencia, RMS y fase en una ventana
común guardada. El ajuste de amplitud de una envolvente que decae se conserva
como diagnóstico adicional; el criterio de amplitud usa la historia y su RMS.

`run_cad_static_convergence(bundlePaths,output,'horseshoe')` y la variante
`'uvlm'` refinan por separado ambos equilibrios. Sus medias no se intercambian:
la dinámica se inicia en su propio equilibrio UVLM. Los anillos sobre la
superficie curva y las herraduras rectas son idealizaciones distintas de la
vorticidad; su diferencia en el CAD es una discrepancia de modelo que debe
informarse, aunque ambas reproduzcan el límite plano del benchmark.

`run_cad_refined_load_check` compara FEM completo y reducción bajo la carga
rígida y la carga en equilibrio deformado. `run_cad_final_review(manifestPath)`
revisa rutas explícitas a las campañas terminadas, guarda tablas, hashes de
evidencia y criterios; no elige automáticamente la carpeta más reciente.
El manifiesto local de la entrega está en `results/cad_original/final_inputs.json`.


## Refinación temporal exacta de operadores

`run_cad_time_refinement(tangentPath)` conserva el equilibrio, CAD y base del
padre y divide dt por dos. `uvlm.refine_time_model` evalúa directamente los
kernels impares y recupera los pares mediante identidades geométricas:

```
W_grueso(a) = W_fino(2a-1) + W_fino(2a)
G_grueso(a) = G_fino(2a) + (G_fino(2a-1) + G_fino(2a+1))/2
```

La segunda identidad mueve los vértices intermedios de cada segmento a la
mitad del movimiento del vértice grueso. También se corrigen D y Lq de la
primera fila de estela; omitir esa corrección cambiaría el problema. Cada
circulación y posición material **fina** conserva su historia independiente:
no se interpola la trayectoria ni se agrupa circulación para ahorrar memoria.

La comparación independiente contra ensamblaje directo en una superficie
curvada da cambios de operadores del orden de 1e−10 o menores y error de
trayectoria 8,8e−13 en el control registrado. El test público contrasta ambos
procedimientos. Los archivos `parent.json` registran hashes del operador y
respuesta originales. Las refinaciones temporales de velocidades distintas se
evalúan por separado de la comparación espacial a dt fijo.

## Presión y movimiento sincronizados

`fsi.dynamic(...,struct('storePressureHistory',true))` conserva la presión
normal equivalente por panel y paso. Un presupuesto explícito de bytes evita
reservas ilimitadas. Los controles que añaden masa inicial pueden proporcionar
`pressureAccelerationLoad` para observar también ese término Bernoulli inicial.

`run_cad_pressure_history(bundlePath)` ejecuta el control no lineal CAD 8×6 y
guarda presión, TE, energía/trabajo/disipación y coordenadas/áreas de paneles
con el mismo tiempo. Es un control pequeño del intercambio, no Cp convergido
de ambas pieles. La variante fina principal conserva su alcance de superficie
equivalente y su propio estudio de discretización.
