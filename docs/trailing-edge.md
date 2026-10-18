# Deformación y vibración del borde de salida

Esta ruta añade una placa delgada y una interfaz estructural común para estudiar
la cola antes del acoplamiento fluido–estructura. Conserva el FEM sólido del CAD.
La placa es una idealización plana: no transforma automáticamente las pieles,
nervaduras y slider del STEP en una estructura equivalente.

## Ejecutar

Desde la raíz del proyecto, con MATLAB base R2026a:

```matlab
addpath('matlab')
run_verification();
comparison = run_tail_benchmarks();
result = run_tail_response();
```

El primer ejemplo compara nueve geometrías publicadas sin ajustar parámetros.
El segundo genera tres figuras MATLAB: desplazamiento estático de la cola a lo
largo de la envergadura, historia temporal central y amplitud frente a frecuencia.
CSV, PNG y MAT se guardan en `results/tail_benchmarks/` y `results/tail_response/`.
Repetir la ejecución sobrescribe las salidas.

## Formulación de placa propia

`plate.mitc4` implementa una placa plana isotrópica de Reissner–Mindlin con
interpolación mixta de cortante, siguiendo a
[Bathe y Dvorkin, 1985](https://doi.org/10.1002/nme.1620210213).
Los grados de libertad son `[w1,theta_x1,theta_y1,w2,theta_x2,theta_y2,...]`.
Theta son pendientes generalizadas: el cortante es
`[dw/dx-theta_x; dw/dy-theta_y]`. Son distintos de los desplazamientos
cartesianos del sólido y de los bloques de coordenadas del modelo `pae`.

Las curvaturas son `[theta_x,x;theta_y,y;theta_x,y+theta_y,x]`. Las matrices
constitutivas son `Db=E*t^3/(12*(1-nu^2))*[1 nu 0;nu 1 0;0 0 (1-nu)/2]` y
`Ds=(5/6)*G*t*I`. El cortante covariante xi se interpola entre los centros de
las aristas inferior/superior, y eta entre izquierda/derecha. El jacobiano
transforma estas componentes al punto de integración. Rigidez y masa se
integran en 2×2 puntos de Gauss. La masa consistente incluye `rho*t`
transversal e inercia rotatoria `rho*t^3/12`. No se utiliza integración de un
punto con modos espurios libres.

El espesor es constante por elemento y la malla debe resolver sus saltos.
`ael.build_plate_case` introduce una línea en cada escalón. Las coordenadas
aumentan en x/y y los nodos del cuadrilátero siguen orden antihorario.
Longitudes en metros, E en Pa y densidad en kg/m³.

La formulación es lineal, sin membrana, curvatura inicial, anisotropía de placa,
pretensión, grandes rotaciones, contacto o fricción. Una piel curva o un
elastómero requiere otra formulación o un equivalente identificado. La placa
evita el bloqueo severo de Tet4 en los benchmarks de flexión delgada; no
resuelve por sí sola la convergencia del STEP.

## Presión distribuida y sensores

`plate.assemble` devuelve `system.P`, que convierte una presión constante por
elemento en fuerzas nodales consistentes. `plate.pressure_map(mesh,xEdges,yEdges)`
acepta otra partición rectangular del campo de presión. Integra la intersección
de cada parche con cada elemento mediante funciones de forma y Gauss 2×2.
Conserva fuerza transversal y primeros momentos de ese campo aun al refinar
la estructura. La huella tiene tamaño físico, independiente del número de vecinos.

Las columnas siguen x más rápido que y, como los paneles VLM. Una presión
positiva es `p_inferior-p_superior` y actúa en +w. Los parches deben estar
dentro de la placa. No se admiten superficies curvas o parches arbitrarios.
`plate.probe` interpola desplazamiento o pendiente en puntos físicos fijos.

El ejemplo aplica la presión equivalente VLM como uniforme sobre cada panel.
Conserva la fuerza vertical, pero el momento cambia: VLM aplica fuerzas en
cuarto de cuerda y esta presión actúa en el centro del panel. Ambos momentos
figuran en `summary.csv`. No se afirma conservar el momento original ni obtener
Cp de las dos caras. La transferencia aerodinámica necesita verificación antes
del acoplamiento cuantitativo.

## Interfaz común para sólidos y placas

La interfaz `ael` recibe matrices y operadores explícitos:

`M*u'' + C*u' + K*u = B*p(t)`, con `y=O*u`.

B transforma amplitudes de presión/fuerza en cargas generalizadas y O
transforma grados de libertad en medidas. Se puede sustituir la placa por
un sólido, con sus propios apoyos y operadores, sin reinterpretar sus DOFs.

```matlab
% s: plate.assemble; fixed: grados de libertad empotrados.
O = plate.probe(s.mesh,[xTE,ySensor],1);
model = ael.reduce(s.K,s.M,fixed,24,0.02,s.P,O);
static = ael.static(s.K,fixed,s.P*pMean,O);
frf = ael.harmonic(model,frequenciesHz,pAmplitude);
```

`ael.reduce` resuelve el problema generalizado de autovalores con apoyos
homogéneos, normaliza en masa y proyecta K, M, B y O. C modal usa razones de
amortiguamiento proporcionadas; el 2% del ejemplo es supuesto. Se guardan
residuos en `eigenResidual`. Comprobar número de modos y flexibilidad estática
antes de interpretar una respuesta reducida.

Para un sólido, construir B con tracciones consistentes o una transferencia
verificada. Si H interpola desplazamientos sólidos en sensores, una matriz O
con 3N columnas y `O(:,3:3:end)=H` mide la componente z. No pasar directamente
esta base intercalada a `pae`, cuyo orden de coordenadas es diferente.

## Presión temporal y respuesta forzada

Una presión estacionaria produce deformación estática. Las vibraciones
requieren carga temporal, perturbación inicial o un fluido que genere la
variación. La cadena desacoplada admite presiones medidas o procedentes de CFD:

```matlab
% pressure: Ncanales × Ntiempos, Pa, con reloj sincronizado.
inputFcn = ael.pressure_history(tPressure,pressure);
dynamic = ael.transient(model,time,inputFcn,q0,v0);
metrics = ael.response_metrics(dynamic.time,dynamic.Y);
comparison = ael.compare_series(dynamic.time,dynamic.Y,tMeasured,yMeasured,chord);
```

`pressure_history` interpola linealmente y rechaza extrapolación temporal.
`transient` usa Newmark de aceleración media, beta=1/4 y gamma=1/2. Los estados
iniciales son coordenadas modales. Para equilibrio inicial de presión media:
`q0=model.K\(model.B*pMean)` y velocidad cero. Un transitorio de arranque no
demuestra inestabilidad aeroelástica.

`harmonic` usa `Re(Y*exp(i*2*pi*f*t))`: amplitudes pico complejas y fase en
radianes. `response_metrics` entrega media, RMS de fluctuación, pico a pico y
bin dominante FFT con ventana Hann, implementada en MATLAB base. El usuario
selecciona la ventana y elimina el transitorio cuando corresponda. Se exige
muestreo uniforme y se informa resolución frecuencial. La frecuencia dominante
de una vibración forzada no identifica por sí sola un modo natural o flutter.

`compare_series` interpola a los instantes medidos, exige cobertura completa
y entrega sesgo, RMSE, error máximo y RMSE/escala. La escala es explícita,
por ejemplo cuerda para error de desplazamiento/cuerda. No desplaza señales
automáticamente para mejorar fase ni ajusta el material.

## Ejemplo del álabe

Se usan cuerda 247 mm, envergadura 360 mm, zona anterior inmóvil hasta 30%
de cuerda y placa posterior de 2 mm. E=2,11 GPa, nu=0,3 y rho=1180 kg/m³
son parámetros del modelo nominal PLA publicado, no de la pieza real.
VLM estacionario se evalúa a 10 m/s y 3° sobre geometría inicial.

Se añade **artificialmente** una oscilación sinusoidal del 10% de la presión
media a 0,8 veces la primera frecuencia estructural. Verifica presión →
vibración → medidas. No predice fluctuaciones del flujo ni el slider real.

Se contrastan 12/24 modos en la frecuencia impuesta, estática modal frente a
completa y dos pasos temporales; las diferencias están en `summary.csv`.
Falta refinación espacial del propio caso del álabe y comprobar la banda modal
completa para otros espectros. Convergencia en otro caso no sustituye esto.

## Estela y siguiente etapa

Las herraduras VLM actuales son rectas y estacionarias. No calculan convección
temporal de lámina vorticosa, roll-up, difusión, separación o calle de vórtices.
Esta versión no genera figuras de estela interpretables como resultados físicos.

1. Identificar estructura, apoyos, slider y amortiguamiento; resolver
   convergencia del CAD o identificar un equivalente estructural válido.
2. Contrastar presión medida/CFD sobre geometría inmóvil con deformación,
   RMS, fase y espectro de la cola en el modelo desacoplado.
3. Verificar acoplamiento estático en dos vías: geometría actualizada,
   conservación de cargas/trabajo, tolerancias y relajación explícitas.
4. Añadir aerodinámica no estacionaria. UVLM representa memoria circulatoria
   y estela vorticosa de flujo adherido; necesita Kutta, conservación de
   circulación, presión no estacionaria, emisión/convección y refinación temporal.
5. Para estela separada, déficit viscoso y vibración inducida por separación,
   contrastar con CFD viscoso y PIV/anemometría. UVLM no los reproduce por sí solo.

Las funciones nuevas no activan realimentación de presión ni confirman flutter.
El alcance está guardado también en `scope.json` del ejemplo.
