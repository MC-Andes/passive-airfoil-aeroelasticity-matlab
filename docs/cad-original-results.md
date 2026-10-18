# Resultados del CAD original y alcance de aceptación

> Esta ruta corresponde al STEP anterior y a sus hipótesis de PLA/montaje.
> El [estudio actual PETG con servo](morphing-wing.md) usa otro STEP y apoyo;
> conserva esta evidencia como histórica y las capacidades como complementos.


El objetivo se conserva: estructura posterior y vibración del STEP original,
con solver FEM/fluido MATLAB propio. No se sustituyó el mecanismo por una placa.
El usuario autorizó hipótesis geométricas y bibliográficas porque no existen
entradas experimentales confirmadas. Las predicciones son **condicionales** a
[esas hipótesis](../cases/cad_original/assumptions.json); no se presentan como
validación del espécimen.

El [primer bloque](cad-original-history.md) se conserva como diagnóstico
histórico. Sus fallos de convergencia motivaron un remallado del cap original,
no un ajuste del material para coincidir con una respuesta deseada.

## Geometría y transferencia

El STEP mantiene SHA256
`5bacb470df1b1239a0d8af2181f10ca14077a4c0eba3b3531126c8981ea4793e`.
La auditoría encuentra un sólido conectado de 587 caras. Sus dos caps son
planos y las normales de las 585 caras laterales no tienen componente de
extrusión en los puntos controlados. Se malla la sección original completa,
se conservan sus curvas y se extruye con una partición tetraédrica conforme.
No se borraron uniones, huecos o caras pequeñas.

El atlas registra firmas de caras, intersecciones, los elementos más esbeltos
y espesores normales hasta la primera salida del material. En las pieles
posteriores candidatas, las medianas son aproximadamente 1,003 mm arriba y
1,176 mm abajo; los mínimos muestreados rondan 0,554 mm. Los máximos incluyen
encuentros con nervaduras y no se interpretan como un espesor uniforme.

El parche de verificación conserva x/c=0,78..0,84 y eta=0,25..0,75. Su área
proyectada es 0,002667558534 m²: a 100 Pa produce Fz=−0,2667558534 N.
Los nueve sensores materiales usan x/c=0,999 y eta=0,05..0,95. La carga se
integra sobre la piel exterior, sin seleccionar cavidades por proximidad.
La interpolación y su adjunto conservan fuerza, momento y trabajo.

La inversión Tri6 y el recorte curvo se contrastan con resultantes analíticas
de un paraboloide. La transformación X,Z,Y se trata como reflexión; se
recalculan normales exteriores y se verifican momentos axiales.

## Tres niveles estructurales aceptados para la hipótesis A

Fuente: `results/cad_original/campaign_20261017T190836945`, contrato
[extruded-fine.json](../cases/cad_original/extruded-fine.json).

| Nivel | DOFs Tet10 | Tetraedros | Máximo TE a 100 Pa [mm] | Error de volumen CAD |
| --- | ---: | ---: | ---: | ---: |
| 1 | 422 994 | 82 782 | 0,219615 | 0,0699% |
| 2 | 752 268 | 149 580 | 0,223638 | 0,0392% |
| 3 | 1 134 945 | 232 452 | 0,224243 | 0,0228% |

Entre los últimos niveles cambia el vector TE 0,2697% y la energía 0,5262%.
La primera frecuencia cambia 0,1544%; los cinco primeros cambios son menores
de 0,52%. El sexto cambia 1,479%, y se conserva como no convergido al 1%.
La banda declarada para esta campaña es 0–60 Hz; contiene el modo fundamental
20,59097 Hz. El criterio de banda no implica convergencia de todos los modos
por encima de ella.

Los residuos estáticos en norma de energía son del orden de 1e−11, y los
modales, de 1e−8 o menores. Los residuos directos de fuerza/eigenvector,
que alcanzan aproximadamente 1,8e−6/5,0e−6 en el nivel final, se conservan.
Se declaró antes de esta campaña un estimador dual de energía para la
aceptación algebraica, debido a los detalles muy esbeltos del CAD. Es un
estimador a posteriori, no una cota certificada. La
[justificación y límites](numerical-residuals.md) permiten evaluar esa decisión;
no se ocultan los residuos directos. Los balances globales de fuerza/momento
del operador estable están cerca de la precisión de máquina.

La proyección de nodos intermedios a curvas/caras del CAD no invierte elementos
controlados. El volumen integrado curvo difiere del BREP aproximadamente
0,00015%. Una **solución FEM completa** sobre esa geometría cambia el vector
TE 0,00190% a igual presión y la primera frecuencia 0,00620%. La primera
frecuencia curva es 20,59225 Hz. Su residuo estático de energía es 5,81e−7.
Un ensayo Ritz inicial con la base geométrica congelada sobreestimó la rigidez;
se conserva como diagnóstico descartado y no como cota del error geométrico.

La cuadratura Duffy 4/6 se comparó en los 232 452 elementos de la geometría
curva, con sus seis modos completos y su desplazamiento estático. El mayor
cambio de energía de esos campos es 8,27e−10 relativo; el cambio del operador
de rigidez proyectado y escalado es 2,87e−10, y el de masa menor de 1e−14.
Pasa la tolerancia de 0,1% para este control. Es convergencia de integración
en esos campos, no una cota de todo el operador ni una segunda solución
completa a orden 6. Se guardan matrices y diagnóstico en
`results/cad_original/meshes/curved_projection/quadrature_comparison.json`.

## Mecanismo y sensibilidades condicionales

La barra interna limitada por las caras 7/9 es la candidata geométrica de guía.
Sus seis coordenadas rígidas medias se obtienen por integración ponderada por
área. La conectividad del sólido se conserva en A/B/C.

| Hipótesis | Restricción adicional | Máximo TE a 100 Pa [mm] | Cambio entre dos últimas mallas |
| --- | --- | ---: | ---: |
| A | Ninguna; sólido monolítico | 0,22424 | 0,270% |
| B | Seis coordenadas de barra bloqueadas | 0,01234 | 1,03% |
| C | Traslación media en cuerda libre; otras cinco bloqueadas | 0,06794 | 0,367% |

Las curvas de carga/descarga elásticas coinciden. Los ciclos adicionales de
Coulomb usan mu=0/0,2/0,4, precarga normal de 1 N y actuación tangencial
hipotética de ±1 N; su disipación es no negativa. Son respuestas de un modelo
de guía supuesto, no histéresis medida. Se guardan también desplazamientos de
guía ±0,1 mm y rigideces finitas de guía, con sus unidades y alcance.

En la malla de sensibilidad, desplazar el límite del empotramiento de 30% a
25%/35% modifica poco la flexibilidad de este cuerpo anterior rígido.
La variante geométrica explícita de altura 0,95/1,05 cambia el máximo TE de
0,21961 mm a 0,24917/0,19464 mm. Esta transformación cambia pieles, nervaduras
y huecos en altura: no es una tolerancia aislada del espesor de pared. El CAD
principal permanece intacto. Las variaciones se separan de las de E, rho,
amortiguamiento y guía.

## Modos y reducción

Las seis frecuencias del nivel final son aproximadamente
20,591; 72,964; 93,597; 98,296; 148,866 y 172,221 Hz. En el fundamental,
81,1% de la energía cinética modal corresponde a elementos con centroide
posterior x/c≥0,70; 87,4% corresponde a desplazamiento en altura. Las otras
componentes y regiones quedan tabuladas, sin clasificar modos por su número.

La base utiliza seis modos FEM y un vector de flexibilidad residual. Para la
carga VLM inicial, el error estático en sensores baja de 0,415% a cerca del
redondeo. La comprobación temporal sobre el **FEM completo** de 1,13 millones
de DOFs usa cuatro períodos y 400 pasos: error máximo en sensores
5,46e−7 respecto a la perturbación inicial y defecto energético 2,07e−7.
Es una comparación en el subespacio del primer modo; el truncamiento bajo
otras cargas se comprueba por separado.

La deformación de una vía a U=10 m/s tiene u/Lf≈0,00149 y rotación máxima
centroidal≈0,00328 rad, por debajo de los disparadores de comparación no lineal
declarados. La formulación no lineal se verifica mediante rotación finita,
tangente consistente y barra de extensión finita; no se afirma una ley
constitutiva calibrada de PLA.

## Aerodinámica, UVLM y controles independientes

La superficie media se deriva de ambas pieles del CAD y mantiene operadores
materiales. El VLM aplica sus fuerzas en el brazo de cuarto de cuerda y se
transfiere por trabajo virtual. A U=10 m/s y alpha=3°, 32×18 paneles dan
CL rígido≈0,21384 y equilibrio de dos vías CL≈0,20963, con máximo TE≈0,21615 mm.
Entre 24×14 y 32×18 paneles cambia aproximadamente 0,2% la deformación;
las mallas más gruesas mostraron diferencias mayores y se conservan.

Theodorsen y Sears verifican heave/ráfaga con k=omega*c/(2U), incluidas
convenciones de fase. En k=0,5/1,6/6,4, las refinaciones finales cumplen
2% de amplitud y 2° de fase. Se estudian por separado paneles, dt, núcleo y
horizonte. La estela conserva todos los anillos emitidos y rechaza exceder
su presupuesto; no pierde circulación por truncamiento.

En movimiento prescrito de la cola real, el UVLM no lineal y su tangente
difieren 5,81e−5 en fuerza perturbada normalizada. En FSI fuerte del CAD,
la diferencia de sensores es 0,0241%; cambiar relajación/tolerancia produce
0,000315%. El trabajo espurio de intercambio es del orden de 1e−25 J y el
defecto energético 3,59e−8. Estos controles aíslan errores del intercambio;
no son validación experimental del amortiguamiento.

La ruta de cálculo y las hipótesis están en [CAD acoplado](cad-coupled.md).
La aceptación final dinámica y viscosa se registra con sus campañas y bandas
en los apartados de cierre que acompañan este informe.


## Rama viscosa: resultados que sí calcula el modelo

La corrección presión/penalización acoplada evita que una proyección posterior
reintroduzca velocidad dentro del cuerpo. Los controles analíticos contrastan
masa, cantidad de movimiento, Poiseuille, Taylor–Green 2D y Beltrami ABC 3D.
En el cilindro Re=20, dominio 48D y malla 768×768, se obtiene Cd=2,04749
frente a 2,05 de Dennis–Chang recogido en la tabla 2 de Russell–Wang.
La diferencia es −0,122%; no se calibró el solver para ese valor ni se infiere
exactitud universal de esa coincidencia. Las corridas anteriores con proyección
separada se mantienen como diagnósticos.

Para el CAD, las cavidades se rellenan solo en la representación **fluida**
del contorno exterior. La variante 3D tapa los extremos; la rama estructural
conserva las cavidades originales. Ambas hipótesis se guardan y no se mezclan.
Los casos 2D positivos usan Re=100, alpha=3°, ventana común t*=12,4–14:

| Malla | Cd | Cl | Cambio del perfil común de déficit |
| --- | ---: | ---: | ---: |
| 384×192 | 0,443290 | 0,152417 | — |
| 512×256 | 0,442268 | 0,148256 | 0,3424% |
| 640×320 | 0,442628 | 0,154457 | 0,4644% |

Los perfiles se comparan en x/c=1,25/1,5/2/2,5 y altura/c=−0,75..0,75.
El Cd cambia 0,0815% entre las dos mallas mayores, pero Cl cambia 4,18%:
la convergencia de estela no demuestra la misma precisión para sustentación.
El paso temporal dividido por dos mantiene las fuerzas prácticamente iguales;
reducir a la mitad el tiempo de penalización cambia Cd≈0,56% y Cl≈1,35%.
El heave impuesto A/c=0,002, omega*=2 es un control de pared móvil, no una
oscilación autoexcitada. Las presiones de piel se muestrean fuera de la capa
penalizada; el offset 1/1,5/2/3 celdas y la referencia del reservorio aguas arriba
quedan registrados. La curva Cp es numérica, no una presión de túnel validada.

Al dimensionalizar Re=100 en este álabe y aire supuesto se obtiene
U≈0,005982 m/s. La tracción viscosa calculada, extruida con fuerza/momento
conservados, produce un máximo TE FEM≈7,61e−11 m. El error de reducción es
0,836%. Ese resultado no se extrapola a U=10 m/s ni prueba que la carga viscosa
sea despreciable al Reynolds del túnel.

El presupuesto para Re≈167166 a U=10 m/s usa una estimación laminar explícita
de capa límite y diez celdas a través de ella: unos 4,02e10 volúmenes uniformes
y 10,3 TB para 32 arrays de trabajo. El equipo dispone de 24 GiB. Este cálculo
motiva el alcance laminar de baja velocidad y el UVLM adherido; no se declara
DNS/LES/URANS del túnel ni se inventa una comparación PIV.

## Comparaciones con literatura

Las doce [fichas](../data/benchmarks/cad_literature.json) distinguen ensayos,
referencias analíticas y referencias de arquitectura. Morettini aporta una
aproximación inicial de material/amortiguamiento de PLA; no identifica Tough PLA
impreso en esta pieza. Theodorsen/Sears y las soluciones viscosas analíticas
verifican ecuaciones. El cilindro es una comparación externa cuantitativa.
La comparación independiente de nueve placas de Huang conserva errores del
primer modo de +0,68% a −10,87%, pese a cambio de malla ≤0,108%; no se ajustó E
ni se afirmó cumplir el objetivo experimental de 5%.

Taguchi usa otra arquitectura de cola y ensayo casi bidimensional a Re≈1,78e5;
la ganancia de sustentación a incidencias mayores no impone el signo de cambio
de CL en este álabe finito a 3°. Kai/Ochi/Ōtomo orientan mecanismos, acoplamiento
y escalas; sin entradas completas y coincidencia de condiciones no se usan
como curvas de calibración. Las hipótesis autorizadas permiten terminar el
estudio numérico y sus bandas, sin convertir referencias ajenas en mediciones
del espécimen.


## Gravedad y ventilación: condiciones de montaje explícitas

La hipótesis primaria coloca la envergadura vertical. Tres cargas gravitatorias
independientes sobre el FEM completo conservan el peso 7,0029 N. La gravedad
en dirección de envergadura da máximo TE en altura de 2,231 micras y máximo
vectorial de 4,223 micras. Los residuos en energía son ≤1,36e−11. La respuesta
al viento se interpreta como incremento sobre la referencia sin flujo, con
estructura lineal y sin corrección de rigidez por pretensión finita.

La orientación importa: gravedad en cuerda/altura produce máximos TE de
0,4458/0,5360 mm. Son otras hipótesis de montaje; no se suman silenciosamente
al caso primario. La orientación real no está medida. Los estados de las tres
orientaciones se guardan en `gravity_orientation/` de la campaña estructural.

La sección original tiene ocho cavidades abiertas en ambos extremos, con
área de puertos por extremo de aproximadamente 0,002957 m². El caso primario
supone presión interna ambiente. Un modelo de impedancia de ventilación,
con gradientes de volumen del sólido real, estima masa añadida del primer modo
entre 0,46% (fuente uniforme en envergadura) y 1,39% (fuente concentrada central).
La presión proxy máxima para TE=10 micras es 0,044 Pa. La resistencia usa una
aproximación hidráulica/oscilatoria, sin medición ni cota certificada.

Si se cerraran los puertos, el cociente proxy de rigidez neumática/estructural
del primer modo sería del orden de 2,96: no puede trasladarse automáticamente
la predicción ventilada a un montaje que selle los extremos. Esa estimación
advierte una condición física concreta; no es una simulación validada de aire
encerrado. La sensibilidad dinámica abierta se conserva por separado de la
rama CFD exterior con extremos tapados.


La recuperación de tracción de piel no tiene la misma precisión que la fuerza
inmersa total: en 640×320 y offset de dos celdas da Cd=0,3646 y Cl=0,1521,
frente a 0,4426/0,1545 por reacción Brinkman. El arrastre de piel queda 17,6%
por debajo; los offsets de una a tres celdas cambian Cd de 0,3769 a 0,3381.
Se conserva ese diagnóstico. La transferencia FEM verifica conservación de
la tracción calculada, sin convertir ese campo en una presión cuantitativamente
validada. Las predicciones de déficit y la reacción global tienen sus propios
controles y no eliminan este error de recuperación local.

El control corto de convección libre del CAD compara Euler/Heun con la estela
congelada. El CL medio Heun cambia 1,62% al dividir dt, y Kelvin se conserva,
pero la geometría de enrollamiento de esa malla gruesa no está convergida.
No se utiliza para la aceptación dinámica: el modelo principal declara
convección ambiental congelada y tiene sus propias refinaciones de memoria,
núcleo, paneles y dt. La circulación calculada no se presenta como déficit viscoso.


La comprobación de formas modales usa cinco estaciones de cuerda por nueve
de envergadura, con las tres componentes. Entre los dos niveles mayores los
MAC diagonales son ≥0,999285 (fundamental 0,999999992). La observación solo en
TE confundía parcialmente dos modos; esa observación limitada se conserva
como diagnóstico y no se usa sola para identificarlos.

En 3D viscoso, las mallas 128×64×64 y 192×96×96 comparten Re=20, dominio y
ventana t*=2,6–3. El perfil medio en el plano central cambia 0,3541%. Cd pasa
de 1,39479 a 1,39021, y Cl de 0,09460 a 0,09965: la diferencia de sustentación
es aproximadamente 5,1% respecto a fina. Son estadísticas de un transitorio
finito, no evidencia de independencia del dominio o de estado estacionario.
El plano transversal x/c≈1,26 muestra déficit, flujo transversal y vorticidad
en corriente calculados, con sus coordenadas y escalas en `cross_wake.mat`.


## Refinación estática ampliada y equilibrio de la dinámica

La secuencia inicial de paneles no bastó. Se añadieron 96, 128, 192, 256,
384 y 512 en cuerda, manteniendo 26 en envergadura y la misma base del sólido.
La geometría exterior tiene un cambio local cerca de x/c≈0,55 que concentra
carga; el muestreo uniforme produce una secuencia no monótona. Las refinaciones
adicionales permiten observar su reducción, sin inferir exactitud de dos
valores coincidentes.

| Modelo estacionario | Paneles en cuerda | Máximo TE reducido [mm] | Cambio vectorial |
| --- | ---: | ---: | ---: |
| Herraduras | 256 | 0,205880 | — |
| Herraduras | 384 | 0,206919 | 0,5024% |
| Herraduras | 512 | 0,207778 | 0,4134% |
| Anillos UVLM | 256 | 0,167449 | — |
| Anillos UVLM | 384 | 0,168577 | 0,6692% |
| Anillos UVLM | 512 | 0,169490 | 0,5385% |

A 512×26, las herraduras dan CL rígido=0,215820 y acoplado=0,211745;
los anillos dan 0,197449/0,194639. El equilibrio UVLM se utiliza para iniciar
la dinámica UVLM. La diferencia entre estas idealizaciones en el CAD es del
orden de 8,1% en CL acoplado y 18,4% en TE. Se informa como discrepancia de
modelo; no se promedian ambos resultados ni se presenta el límite plano del
benchmark como validación de estas cargas sobre el CAD curvo.

Las cargas de herraduras a 512×26 se comprobaron sobre el FEM completo de
1,13 millones de DOFs: máximos TE de una/dos vías 0,211025/0,207174 mm,
con errores de reducción 0,2872%/0,2916% y residuos de energía ≤1,42e−11.
En la segunda comparación se resuelve el sólido completo bajo la carga del
equilibrio deformado; se comprueba la reducción sin reemplazarla por una placa.

El mapa de presión del equilibrio UVLM y el término Bernoulli inicial se
guardan por panel, con áreas y unidades Pa. Este último se obtiene de la
aceleración calculada por el acoplamiento: su proyección modal cierra con
error relativo 8,23e−16. Es una contribución inicial de presión fluctuante,
no la fluctuación total ni Cp independiente de las dos pieles. Los picos
locales próximos al cambio geométrico no se tratan como presiones de túnel
convergidas.


La comprobación FEM completa de cargas de **anillos UVLM 512×26** da máximos TE
rígido/equilibrado de 0,170844/0,168981 mm. El error de la reducción es
0,2970%/0,3011%, con residuos de energía ≤1,56e−11. Esto comprueba ambas cargas
sobre el sólido original, manteniendo separada la discrepancia respecto a
herraduras.

El control no lineal de presión sincronizada usa 8×6 paneles, dt=0,0005 s y
ventana 0–0,1 s. La perturbación máxima de presión equivalente calculada es
0,0369 Pa; el defecto energético es 3,38e−8 y el del trabajo de interfaz
1,55e−25 J. Presión de cada panel, movimiento TE y energía se guardan con el
mismo reloj en `results/cad_original/pressure_history`. Es un control de
acoplamiento en una malla gruesa; no se presenta su distribución como Cp de
túnel ni como solución espacialmente convergida.


## Dinámica final, mapas y bandas

La [ruta acoplada](cad-coupled.md) describe los operadores y sus verificaciones.
La revisión local `results/cad_original/final_review/acceptance.json` acepta
los criterios numéricos declarados y conserva el manifiesto y hashes de evidencia.
`report.md` en esa carpeta integra las figuras calculadas.

En A, la perturbación inicial es 10 micras, sin presión sinusoidal impuesta.
Se usa zeta=0,008 supuesto; la ventana de respuesta es 0–1,2 s. El ajuste excluye
dos picos iniciales y amplitudes menores de 0,1% del máximo. Las frecuencias
se estiman por picos del primer modo, no por el máximo de una FFT mezclada.

| U [m/s] | f acoplada [Hz] | Decaimiento [1/s] | RMS TE máximo [micras] | Picos válidos |
| --- | ---: | ---: | ---: | ---: |
| 2,5 | 20,31255 | −2,43419 | 2,9240 | 22 |
| 5 | 20,33591 | −3,80476 | 2,3444 | 22 |
| 7,5 | 20,37599 | −5,16383 | 2,0155 | 22 |
| 10 | 20,42873 | −6,58110 | 1,7887 | 20 |

Las tres primeras usan dt=31,25 microsegundos; U=10 usa 15,625 microsegundos.
Los defectos de energía son ≤2,63e−14 y los residuos relativos de equilibrio
≤1,59e−13. El RMS corresponde al transitorio completo declarado; no es una
amplitud estacionaria ni una desviación experimental.

En U=10, 96→128 paneles en cuerda cambia frecuencia 0,0222%, decaimiento
1,63%, máximo de historia/pico 0,62%, RMS 0,80% y fase 0,38°.
La división final de dt (31,25→15,625 microsegundos) cambia decaimiento 1,079%,
máximo de historia/pico 0,439%, RMS 0,537% y fase 0,411°. La reducción anterior
62,5→31,25 no cumplía 2% y motivó el nivel adicional; se conserva ese resultado.
Los controles de núcleo 1e−5→1e−6 m y longitud inicial 1000→500 cuerdas pasan
por separado. La comparación 32×18/32×26 controla la discretización en envergadura.

La convergencia temporal también se comprobó en U=2,5/5/7,5: cambios de
decaimiento 0,288%/0,919%/1,608% y fases 0,222°/0,452°/0,648°. Las comparaciones
espaciales de esas velocidades usan dt idéntico; no mezclan error temporal y
espacial. Cada pareja y ventana se conserva en la revisión.

| Guía hipotética a U=10 | f [Hz] | Decaimiento [1/s] | Cambio temporal del decaimiento | Fase |
| --- | ---: | ---: | ---: | ---: |
| B, seis coordenadas bloqueadas | 57,32604 | −4,53228 | 0,262% | 1,649° |
| C, traslación en cuerda libre | 37,93598 | −7,43506 | 0,812% | 0,717° |

B usa una ventana de 0,4 s y C de 0,6 s, ambas con 20 picos válidos en la
corrida final. C necesitó refinación adicional: 125→62,5 microsegundos no
cumplía el criterio; las parejas finales se guardan explícitamente.
No se extrapolan estos resultados a holgura o fricción reales sin identificar.

En la base A, 1/3/6/7 coordenadas dan errores máximos de sensor normalizados
por la perturbación inicial de 0,1272%/0,0275%/0,0213%/0%. La corrección estática
residual y los modos adicionales apenas modifican esta respuesta en la banda
0–60 Hz, aunque siguen siendo necesarios para representar otras distribuciones
de carga y otras bandas.

Las nueve variantes E/densidad/zeta producen máximos TE estáticos UVLM de
0,14365–0,21368 mm y decaimientos de −7,8781 a −5,5889 1/s, con al menos 17
picos válidos. Son una envolvente determinista, sin probabilidad o nivel de
confianza asignados. Se usó dt=31,25 microsegundos y la tangente fluida base;
la refinación temporal adicional de A se cuantifica aparte. La banda no cubre
por completo la discrepancia de idealización aerodinámica ni identifica
propiedades de la impresión real.

El análisis de ventilación abierto cambia la frecuencia de 20,433 a
20,388/20,295 Hz para factores de impedancia 1/3 y 1; el decaimiento pasa de
−6,5101 a −6,5072/−6,5010 1/s. Estos valores son una sensibilidad del modelo
proxy de cavidades, sin inferir una resistencia de ventilación medida.

## Forma, pendiente y curvatura posteriores

Se siguen 64 estaciones materiales x/c=0,68..0,995 y las mismas nueve eta.
Las formas dinámicas combinan el equilibrio UVLM del propio modelo con sus
perturbaciones a t=0/0,0125/0,025/0,05 s. No se imponen esas formas al solver.
Los datos y figuras están en `results/cad_original/posterior`.

La pendiente y la curvatura firmada de la proyección cuerda/altura se obtienen
con un ajuste cúbico en ventana ±0,025c, de longitud total 12,35 mm. Se evalúan
en x/c=0,705..0,970, evitando la arista geométrica y derivadas de facetas.
El contraste con una parábola analítica da error máximo 2,97e−13 1/m.

Bajo el parche de 100 Pa, el máximo cambio de pendiente ajustada es 0,003174 rad
y el de curvatura 0,08928 1/m en la malla final. El campo de curvatura cambia
3,20% entre las dos últimas mallas; su precisión es inferior a la de TE y
energía. Se conserva como estimación geométrica en esa ventana, sin declarar
curvatura puntual, tensión o resistencia convergidas al 2%.

## Escalas e interpretación física

Re en las cuatro velocidades es aproximadamente 41791/83583/125374/167166;
b/c=1,45751. Con k=pi*f*c/U se obtiene 6,3047/3,1560/2,1081/1,5852. Los
benchmarks no estacionarios abarcan ese intervalo de frecuencias reducidas.
La razón global m_s/(rho_aire*S*c)=26,533 incluye el material del cuerpo anterior
sujeto; no es una razón de masas añadidas del primer modo. La flexibilidad TE
bajo el parche fijo es aproximadamente 0,8406 mm/N, y cambia con la distribución
de carga; no se asigna una rigidez D de placa arbitraria al mecanismo.

Con las hipótesis estudiadas, las respuestas A/B/C son transitorios que decaen.
No se obtuvo una oscilación autosostenida ni se estableció una frontera real de
flutter. La presión calculada, el movimiento y los balances respaldan esa
clasificación dentro de la banda y ventana declaradas. Los datos del espécimen,
Cp de túnel y PIV no se han sustituido por figuras de otra geometría ni por
mediciones inventadas.
