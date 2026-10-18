# Plan de desarrollo aeroelástico sobre el CAD original

> Esta ruta corresponde al STEP anterior y a sus hipótesis de PLA/montaje.
> El [estudio actual PETG con servo](morphing-wing.md) usa otro STEP y apoyo;
> conserva esta evidencia como histórica y las capacidades como complementos.


Fecha de referencia: 17 de octubre de 2026. Base de código inspeccionada:
commit `9d77047876ba6705589eac924fccf1424606448c`, rama
`codex/trailing-edge-benchmarks`, [PR 4](https://github.com/MC-Andes/passive-airfoil-aeroelasticity-matlab/pull/4).
La base inicial se conserva como referencia. El bloque de ejecución del
17 de octubre de 2026, registrado al final, registra el desarrollo numérico completado, sus criterios de aceptación y
los límites de validación física que las hipótesis no convierten en mediciones.

## Objetivo y resultado esperado

Desarrollar en MATLAB un solver propio que estudie el álabe del archivo
**Airfoil with slider alone.stp**, fijo en una corriente de aire. La prioridad
es predecir la forma y vibración de su región posterior, explicar qué cargas
de presión las producen y estudiar la estela con un modelo fluido apropiado.

La geometría original será el caso principal. Los modelos de placa, viga y
perfiles publicados servirán para verificar componentes, identificar límites
y contrastar tendencias. No sustituirán automáticamente el mecanismo del CAD.
El desarrollo seguirá la secuencia estructura → cargas sobre forma fija →
respuesta desacoplada → acoplamiento estático → aerodinámica no estacionaria →
acoplamiento dinámico. La estela viscosa constituye una extensión específica.

El resultado final deberá incluir entradas reproducibles, geometría y regiones
trazables, convergencia, incertidumbre, comparaciones justificadas y figuras
MATLAB sencillas. Coincidir visualmente con una figura de otro perfil no será
evidencia suficiente de validez.

## Estado inicial comprobado

| Elemento | Estado al iniciar el plan |
| --- | --- |
| Archivo fuente | STEP aportado por el usuario; copia local original y referencia disponible en Downloads. El Inventor IPT queda como fuente auxiliar. |
| SHA256 STEP | `5bacb470df1b1239a0d8af2181f10ca14077a4c0eba3b3531126c8981ea4793e` |
| Topología importada | Un sólido BREP válido y 587 caras; esto no define contactos, materiales o grados de libertad del slider. |
| Dimensiones exteriores aproximadas | Cuerda 246,996 mm, envergadura 360 mm y altura máxima de la caja envolvente 30,394 mm. Esta altura no es el espesor de la piel. |
| Ejes CAD | X en cuerda, Z en envergadura, Y en altura. La ruta actual usa componentes canónicas X,Z,Y. |
| Volumen auditado | 0,0006155285293484151 m³, integración adaptativa OpenCASCADE. |
| Malla gruesa Tet4 | 21 149 nodos, 69 082 elementos, error de volumen 0,476%. |
| Malla fina Tet4 | 85 367 nodos, 286 037 elementos, error de volumen 0,102%. |
| Calidad de malla | Existen elementos muy esbeltos; la razón mínima de valores singulares del jacobiano es del orden de 1,6×10⁻⁶. Se debe localizar dónde aparecen. |
| Referencia estructural histórica | Tough PLA isotrópico supuesto, región anterior x/c≤0,30 inmóvil, U=10 m/s, alpha=3°, rho=1,225 kg/m³. No representa todavía un montaje identificado. |
| Desplazamiento máximo histórico | 0,005520 mm y 0,016286 mm en gruesa/fina; cambio relativo a fina de aproximadamente 66,1%. No son necesariamente desplazamientos del borde de salida. |
| Tensión máxima histórica | 0,374 y 1,647 MPa; valores locales dependientes de malla y carga. No utilizarlos como resistencia validada. |
| Solver existente | FEM Tet4 propio, masa consistente, materiales isotrópicos/ortótropos, estática, modos, Newmark y VLM estacionario. |
| Ruta de cola añadida | Placa plana MITC4, presión prescrita, adaptador modal y comparación de nueve placas PLA publicadas. No representa aún la estructura interna del CAD. |
| Verificación | 33 controles numéricos pasan en MATLAB R2026a; los ejemplos nuevos usan MATLAB base. |
| Capacidades pendientes | Convergencia estructural del CAD, mecanismo del slider, grandes deformaciones del sólido, UVLM, estela viscosa y acoplamiento dinámico del CAD. |

Los valores se trazan a `data/reference/cad_info.json`, `cad_volume_audit.json`,
`data/fem_reference/mesh_info.json`, `data/fem_refined/mesh_info.json` y
`run_fem_convergence`. El volumen no adaptativo histórico, aproximadamente
0,0005430 m³, no debe emplearse como volumen exacto: su integración produjo
un error que ya se corrigió en la auditoría.

## Preguntas físicas que debe resolver el estudio

1. ¿Qué parte de la cola se deforma, cuánto y en qué dirección bajo una carga conocida?
2. ¿Qué modos movilizan el borde de salida y cómo cambian con el slider y el montaje?
3. ¿Qué presión media y fluctuante actúa sobre cada cara y sobre qué área?
4. ¿La vibración es un transitorio, respuesta a excitación externa o inestabilidad del sistema acoplado?
5. ¿Se necesita la estela circulatoria de flujo adherido o desprendimiento y déficit de velocidad viscosos?

Una carga estacionaria aplicada a una estructura lineal amortiguada produce
un equilibrio tras el transitorio. No atribuir una oscilación sostenida al
flujo si la presión oscilante se introdujo manualmente.

## Información experimental y escenarios provisionales

Solicitar al equipo experimental: material y proceso exactos, orientación y
relleno, masa y espesores medidos, propiedades direccionales, curvas de carga,
amortiguamiento, montaje y rigidez del soporte, estado inicial, funcionamiento
y bloqueo del slider, velocidades/incidencias y paredes del túnel. Para dinámica,
pedir presiones y deformaciones sincronizadas, posiciones de sensores, ancho de
banda, muestreo, incertidumbre y repeticiones. Para estela, solicitar planos y
resolución PIV o puntos y respuesta de anemometría. La lista extensa está en
[datos experimentales](experimental-data.md).

Mientras llegan esos datos, continuar las tareas geométricas, los benchmarks,
la convergencia y las sensibilidades. Usar escenarios explícitos:

| Escenario | Representación permitida | Interpretación |
| --- | --- | --- |
| A monolítico | Sólido conectado tal como se importó; material y apoyos declarados. | Referencia numérica de esa idealización, no demostración de cómo funciona el slider. |
| B slider bloqueado | Restricciones relativas y posición fija, cuando cuerpos/superficies estén identificados. | Predicción para una configuración mecánica documentada. |
| C slider guiado | Movimiento permitido en su dirección física; topes, fricción y precarga solo cuando se conozcan. | Modelo del mecanismo móvil y sus incertidumbres. |

No separar nodos, inventar bisagras ni eliminar uniones del sólido únicamente
para obtener mayor flexibilidad. Si el archivo no permite representar el
movimiento real, solicitar el ensamblaje o una explicación del mecanismo y
conservar A como hipótesis mientras se prepara la interfaz de restricciones.

## Magnitudes comunes para todos los casos

Definir en coordenadas materiales estaciones de envergadura eta=y/b y
posiciones x/c. Seguir los mismos puntos físicos tras deformación y refinamiento.

- Desplazamiento firmado del borde de salida, `w_TE(eta,t)`, y normalizado `w_TE/c`.
- Forma posterior `w(x/c,eta,t)` y cambio de curvatura; pendiente de cola obtenida
  sobre un intervalo definido, sin diferenciar una arista geométricamente singular.
- Frecuencias, formas modales, participación de la cola y amortiguamiento.
- Media, RMS fluctuante, pico a pico, amplitud/fase y espectro en una ventana declarada.
- Cp superior/inferior cuando el fluido los resuelva; diferencia de presión,
  fuerzas y momentos integrados, CL y Cm con las mismas referencias.
- Para estela: circulación y geometría vorticosa en UVLM; velocidad media,
  déficit, fluctuaciones y vorticidad solo cuando el modelo las represente.

Registrar `Re=rho*U*c/mu`, relación b/c, fracción flexible, razón de masas con
su definición, y `k=omega*c/(2*U)=pi*f*c/U`. Si una fuente define k de otra
forma, convertirlo explícitamente. Strouhal `St=f_s*c/U` requiere definir f_s
y la longitud usada. No trasladar sin conversión una escala hidrodinámica a aire.

Para placas equivalentes puede usarse `Ca=q*Lf^3/D`, con q=rho*U²/2 y D la
rigidez de placa por unidad de ancho. Para el CAD complejo no asignar un D
arbitrario: comparar primero una flexibilidad efectiva `C_TE` medida/calculada
en m/N bajo una distribución normalizada de carga y conservar esa distribución.

## Secuencia de trabajo y dependencias

Las etapas 1 y 2 alimentan la etapa 3. Las etapas 3 a 5 deben cerrar el modelo
estructural. La etapa 6 puede avanzar en paralelo con sus verificaciones, pero
la etapa 7 requiere ambas rutas. Las etapas 8 a 10 dependen del resultado de 7.
La ruta viscosa de 11 responde a una pregunta distinta de la estela UVLM.

### Etapa 1 Congelar y reproducir la referencia

**Trabajo.** Registrar versión de MATLAB, commit, hashes, unidades, material,
apoyos, discretización aerodinámica y transferencia. Ejecutar verificación y
referencia actual. Separar configuración física de opciones numéricas y rutas.
Crear un contrato de caso con `geometry`, `regions`, `material`, `supports`,
`slider`, `flow`, `loads`, `sensors`, `discretization` y `provenance`.
Marcar cada dato como medido, publicado o supuesto.

**Entregable.** Configuración del caso CAD y reporte reproducible de las dos
mallas, conservando la corrida histórica. Resultados nuevos en otra carpeta.

**Cierre.** Las 33 verificaciones siguen pasando y la referencia se reproduce
sin archivos personales salvo el STEP cuando se regenere geometría. Registrar
la dependencia de herramientas de preprocesamiento sin delegarles el FEM.

### Etapa 2 Identificar regiones y superficies del CAD

**Trabajo.** Inspeccionar secciones y caras en múltiples estaciones. Etiquetar
piel superior/inferior, borde de salida, zona anterior, nervaduras, flexuras,
slider, interfaces y apoyos, según evidencia geométrica y experimental.
Medir espesores locales, huecos, radios y longitudes flexibles. Crear un mapa
que vincule entidades CAD y grupos físicos; los números de cara pueden cambiar
al remallar, por lo que se deben acompañar de firmas geométricas y comprobaciones.

Distinguir frontera del sólido y superficie mojada por el aire. Las paredes
de huecos internos no se cargan automáticamente como piel exterior. Si existen
aberturas que comunican con el flujo, decidir explícitamente cómo se modela
la presión interna. Revisar los picos de las envolventes actuales cerca del
mecanismo y de los extremos; cinco secciones son una representación preliminar.

Documentar la transformación CAD–solver. La permutación X,Z,Y tiene determinante
negativo: orientar de nuevo las caras con su volumen propietario, transformar
correctamente momentos y bases materiales, y verificar signos con una carga conocida.

**Entregable.** Atlas geométrico, tabla de espesores, grupos físicos, sensores
TE y transformación verificable. Se conserva el STEP sin modificaciones ocultas.

**Cierre.** No hay caras cargadas o apoyadas por cercanía accidental a otra
pieza; cada región importante tiene función y procedencia. Las ambigüedades
del slider quedan registradas como escenarios, no como conectividad inferida.

### Etapa 3 Corregir la transferencia antes de atribuir todo a Tet4

**Trabajo.** Separar convergencia estructural de aplicación de carga. Primero
aplicar una tracción o diferencia de presión analítica sobre un parche de
tamaño físico fijo, en la piel identificada. Repetir en ambas mallas con el
mismo apoyo. Después estudiar cargas VLM sin cambiar simultáneamente su malla.

La transferencia actual usa 16 vecinos por punto por defecto. Aunque reproduce
trabajo/fuerza/momento, puede contraer su huella al refinar o cruzar huecos.
Implementar integración sobre caras físicas y una transferencia de superficies
que respete regiones. Para una misma interpolación cinemática H, mantener
la relación adjunta de cargas `F_s=H'*F_a`, considerando rotaciones/offsets si
se utilizan láminas. Una presión medida se integra en su área, no se convierte
sin control en una fuerza puntual.

Verificar fuerza, momento respecto a un origen común y trabajo virtual. No
confundir la conservación de una presión uniforme por panel con el momento
de una fuerza VLM aplicada en cuarto de cuerda: actualmente difieren.

**Entregable.** Diagnóstico que cuantifique cuánto cambia la flexibilidad por
malla, huella y región de transferencia. Nuevo operador para pieles curvas.

**Cierre.** La carga y el área física permanecen iguales entre refinaciones.
Las diferencias restantes pueden investigarse como discretización estructural
y geometría, sin atribuirlas prematuramente a una sola causa.

### Etapa 4 Implementar elementos adecuados conservando el sólido

**Ruta principal.** Añadir Tet10 a la geometría original. Primero verificar
tetraedros de aristas rectas; después evaluar nodos intermedios ajustados al
CAD para geometría curva. Elevar el orden sobre una malla Tet4 no equivale a
mejorar automáticamente la representación geométrica o la calidad de elementos.

Implementar funciones de forma, gradientes, B, rigidez, masa consistente y
cuadratura con orden suficiente para cada integral. Para geometría curva,
comprobar convergencia de cuadratura y positividad del jacobiano en puntos
de integración y puntos de control adicionales.

Adaptar ensamblaje, caras Tri6, orientación, tracciones, extracción de frontera,
conectividad, apoyos, recuperación de tensiones y posprocesado. Cambiar solo
el archivo del elemento dejaría inconsistentes el resto de operadores.
Verificar movimientos rígidos, deformación afín, masa, permutación de nodos y
aristas, flexión/torsión y modos contra referencias independientes.

Localizar y corregir elementos esbeltos próximos a flexuras y encuentros.
Usar tamaños ligados al espesor y curvatura, controles locales y varios niveles
de malla. Eliminar una cara pequeña solo si se demuestra que no altera una
unión o flexibilidad relevante; toda simplificación produce una variante trazable.

**Alternativa condicionada.** Si el coste del sólido no permite estudiar las
pieles finas, extraer superficies medias del CAD y desarrollar láminas curvas
con espesor, membrana, flexión, orientación y acoplamiento a sólidos. La placa
MITC4 plana actual no es automáticamente una lámina curva. Validar conexiones,
offsets y energía frente a submodelos sólidos antes de usar el modelo híbrido.

**Entregable.** Tet10 verificado y comparación Tet4/Tet10 sobre las mismas
regiones y cargas. Decisión documentada sobre sólido completo o híbrido CAD.

**Cierre.** Tres niveles de discretización estructural y convergencia de las
medidas de cola, energía y modos relevantes. Mantener tensiones singulares
fuera del criterio principal de desplazamiento.

### Etapa 5 Identificar apoyos material y mecanismo

**Trabajo.** Empezar con ensayos estáticos sin flujo: cargas conocidas en áreas
definidas, varios niveles, posiciones del slider y descarga. Comparar forma
completa además de un único desplazamiento. Incorporar masa y gravedad cuando
la orientación y la deformación inicial del ensayo lo requieran.

Usar propiedades publicadas como rangos iniciales, nunca como caracterización
del filamento real. Introducir ortotropía con ejes de impresión identificados.
Si el relleno se homogeniza, comprobar masa, rigidez y escala de validez.
Separar sensibilidad a espesor, apoyo, E y slider antes de ajustar parámetros.

Incorporar restricciones multipunto o multiplicadores para guías/uniones
documentadas. Añadir contacto, topes y fricción después de verificar sus casos
elementales. Un slider libre puede introducir movimientos rígidos físicos:
no eliminarlos artificialmente con resortes numéricos sin reportarlo.

Para deformaciones apreciables, incorporar cinemática no lineal y tangente
consistente antes de interpretar la respuesta acoplada. Como señales de revisión,
usar desplazamiento/Lf del orden de 5% o rotaciones del orden de 0,1 rad;
son disparadores de comparación lineal–no lineal, no límites universales.
Elastómeros requieren una ley material identificada y posiblemente viscoelasticidad.

**Entregable.** Escenarios A/B/C con propiedades y restricciones explícitas,
curvas carga–deflexión, histéresis cuando exista y sensibilidades.

**Cierre.** El ajuste no compensa una malla no convergida o un apoyo incorrecto.
Si faltan ensayos, presentar bandas condicionales de predicción y mantener
pendiente la validación del mecanismo.

### Etapa 6 Modos y dinámica estructural del CAD

**Trabajo.** Calcular modos desde K/M del modelo convergido y, si corresponde,
desde la tangente alrededor de su equilibrio precargado. Identificar flexión
en cuerda, envergadura, torsión y movimientos locales de piel/slider. Clasificar
por forma y participación de la cola, no por número de modo solamente.

Comparar modos entre mallas en sensores comunes y usar MAC con ponderación
definida. Modos próximos o repetidos requieren comparar subespacios. Medir
amortiguamiento por decaimiento o FRF. Extender el adaptador `ael.reduce`
solo después de verificar los operadores de observación del sólido.

Comparar respuesta completa y reducida bajo cargas conocidas, flexibilidad
estática y convergencia del número de modos. Considerar corrección estática
residual si los modos retenidos no representan adecuadamente las cargas.

**Entregable.** Primeros modos relevantes del CAD, FRF y respuesta libre,
observables TE y modelo reducido con errores cuantificados.

**Cierre.** No hay modos espurios ni pérdida de masa/apoyos en la reducción.
La banda de frecuencia de interés está resuelta en espacio, modos y tiempo.

### Etapa 7 Aerodinámica fija y respuesta desacoplada del CAD

**Trabajo aerodinámico.** Construir la superficie del álabe desde la piel
identificada, conservar su forma real y verificar secciones, extremos y área.
Usar VLM como nivel inicial de cargas integradas en flujo adherido. Refinar
paneles y longitud de herraduras independientemente de la estructura.
Evaluar si la superficie media representa el objetivo de presión solicitado.

El espesor exterior del álabe y el Reynolds del túnel pueden hacer relevantes
los efectos viscosos. Para Cp de ambas caras, separación o presión en cavidades,
se necesita un modelo que los resuelva o datos medidos. Las correcciones
empíricas deben identificarse y validarse; no presentar VLM como CFD viscoso.

**Trabajo desacoplado.** Primero usar presión conocida sobre forma inicial;
después presión temporal medida o de una referencia independiente. Integrar
la carga en las caras CAD y calcular w_TE, giro, media, RMS, fase y espectro.
Guardar presión, movimiento y energía con el mismo reloj.

**Entregable.** Primera campaña CAD desacoplada con cargas físicas fijas,
forma deformada, perfiles TE, modos, FRF y balances. Las excitaciones sintéticas
se mantienen como controles del solver y se etiquetan como impuestas.

**Cierre.** Estructura y cargas tienen verificaciones independientes. La
diferencia con literatura se interpreta considerando geometría, Reynolds,
apoyos y propiedades, sin forzar coincidencia de una única curva.

### Etapa 8 Acoplamiento estático en dos vías

**Trabajo.** Actualizar superficie y normales con la deformación, resolver
presión y transferir cargas de vuelta. Mantener correspondencia material de
sensores y caras. Aplicar continuación en presión dinámica/incidencia,
relajación y límites de iteración; informar fallo cuando no haya convergencia.

Comparar el límite de alta rigidez con el caso rígido y el límite de carga
pequeña con una vía. Verificar transferencia y signos en geometría deformada.
La carga seguidora y la cinemática no lineal deben ser coherentes cuando se
incluyan. Separar residual de iteración y cambio por refinación.

**Entregable.** Equilibrios acoplados del CAD, historia de residuos, diferencia
frente a una vía y comparación forma–carga con referencias compatibles.

**Cierre.** El equilibrio es independiente de tolerancias y relajación dentro
de la precisión exigida; no se interpreta falta de convergencia numérica
como divergencia física sin análisis adicional.

### Etapa 9 Fluido no estacionario y estela circulatoria

**Trabajo.** Implementar UVLM propio en MATLAB sobre forma rígida primero:
anillos ligados, condición de Kutta, emisión y convección temporal de estela,
conservación de circulación y fuerzas con contribución no estacionaria.
No inferir toda la presión dinámica de la fórmula estacionaria de sustentación.

Verificar límite estacionario, movimiento armónico pequeño y ráfagas con
referencias analíticas o numéricas independientes bajo iguales hipótesis.
Estudiar paso temporal, paneles, longitud de estela, núcleo vorticoso y esquema
de convección. Conservar circulación al truncar o agrupar la estela.

Pasar después a movimiento prescrito de la cola CAD, aislando errores
aerodinámicos antes de permitir respuesta estructural libre. Distinguir
filamentos/láminas vorticosas calculadas de líneas de corriente decorativas.

**Entregable.** Cargas no estacionarias verificadas y estela circulatoria con
su estudio de resolución. Se declara el régimen de flujo adherido aplicable.

**Cierre.** Los errores de tiempo, estela y discretización no dominan amplitud
y fase. UVLM no se usa para afirmar déficit viscoso o desprendimiento por separación.

### Etapa 10 Acoplamiento dinámico del CAD

**Trabajo.** Acoplar el fluido verificado y la estructura mediante intercambio
consistente de posición, velocidad y carga. Empezar cerca del equilibrio de 8,
con perturbación pequeña y amortiguamiento identificado. Si la interacción
fluido–estructura lo exige, iterar dentro de cada paso temporal hasta converger;
un esquema débil de un intercambio necesita demostrar estabilidad y precisión.

Comprobar subiteraciones, número de modos, paso temporal, conservación de trabajo
y balance entre energía estructural, trabajo aerodinámico y disipación. Usar
varias velocidades y perturbaciones para separar respuesta forzada, decaimiento,
crecimiento e inestabilidad numérica. Una simulación breve no establece flutter.

**Entregable.** Historias CAD acopladas y estela circulatoria, mapas de respuesta
con bandas de incertidumbre y límites del régimen modelado. Una oscilación
autosostenida de amplitud finita exige no linealidades físicas adecuadas.

**Cierre.** Frecuencia, amplitud, fase y crecimiento/decaimiento son estables
ante refinación y subiteraciones. El mecanismo que excita la vibración queda
identificado y no depende de una presión sinusoidal añadida sin origen físico.

### Etapa 11 Estela viscosa cuando la pregunta experimental lo requiera

Si el objetivo incluye separación, desprendimiento vorticoso, déficit de
velocidad o fluctuaciones turbulentas que excitan la cola, UVLM resulta
insuficiente. Esta rama requiere Navier–Stokes y es la ampliación más costosa.

Manteniendo el desarrollo propio, verificar primero un solver viscoso 2D
MATLAB con conservación de masa, canal analítico, vórtice transitorio y un
caso externo publicado. Estudiar luego secciones extraídas del CAD con pared
móvil y transferencia presión/tracción. Un resultado 2D no representa los
extremos ni los vórtices de punta del álabe completo.

La extensión al CAD 3D exige estrategia de malla o frontera inmersa, movimiento
de pared, presión, condiciones del túnel, tratamiento laminar/transicional o
turbulento y convergencia espacial/temporal. No proponer DNS a Reynolds de túnel
sin estimar resolución, memoria y coste. URANS/LES requieren modelos y
validación propios; no son una opción que se activa sin desarrollo adicional.

Los solvers externos pueden aportar referencias independientes, sin sustituir
silenciosamente el núcleo propio solicitado. Si se decide usarlos como motor
de resultados, será un cambio explícito de alcance del proyecto.

**Entregable.** Predicciones de estela correspondientes al modelo viscoso
efectivamente implementado y comparación PIV/anemometría en los mismos planos.
Esta rama puede avanzar después del primer resultado estructural y UVLM útil.

## Criterios cuantitativos propuestos

Son objetivos iniciales de ingeniería, no tolerancias tomadas de la literatura
ni garantía universal. Ajustarlos antes de las campañas según magnitud,
incertidumbre experimental y coste; registrar cualquier cambio.

| Comprobación | Criterio inicial de cierre |
| --- | --- |
| Geometría y masa | Volumen a menos de 1% del CAD auditado; masa comparada con medición cuando exista y tolerancias físicas declaradas. |
| Equilibrio y operadores | Residuos algebraicos normalizados del orden de 10⁻⁶ o menores, con escalas físicas que eviten dividir por una carga casi nula. |
| Transferencia | Defectos relativos de fuerza, momento y trabajo de hasta 10⁻⁶ en pruebas controladas; conservar la misma área y región al refinar. |
| Estructura CAD | Cambios menores de 2% en w_TE/energía y menores de 1% en frecuencias relevantes entre los últimos niveles; comprobar al menos tres niveles y tendencia. |
| Tensiones | Evaluar valores en regiones y distancias físicas definidas; objetivo 5% si esa magnitud es necesaria. No exigir convergencia del máximo en una singularidad. |
| Reducción y tiempo | Cambios menores de 2% en amplitud/RMS y 2° en fase en la banda objetivo; primera estimación dt≤1/(50*f_max_interes), seguida de refinación real. |
| Acoplamiento estático | Residuos relativos de carga y desplazamiento de hasta 10⁻⁴, más independencia de tolerancia; declarar escalas de carga y desplazamiento. |
| Acoplamiento dinámico | Cambios menores de 2% al refinar tiempo/subiteración y balance energético comprobado; cuantificar trabajo espurio del intercambio. |
| Comparación experimental | Objetivo inicial de 5% en frecuencias y 5–10% en deformación cuando los datos lo permitan; evaluar incertidumbre combinada. No declarar cumplimiento si el error excede el presupuesto. |

Una diferencia pequeña entre dos mallas no demuestra por sí sola exactitud.
Si la secuencia no es monótona o está dominada por mala geometría/contacto,
añadir niveles o cambiar discretización. Los errores numéricos deben quedar
por debajo del efecto físico que se pretende distinguir.

## Literatura y función de cada referencia

| Referencia | Uso concreto | Condición para una comparación válida |
| --- | --- | --- |
| [Huang y Zhang 2026, datos v2](https://data.mendeley.com/datasets/8jxw5533py/2) | Verificar flexión/torsión y transición rígida–flexible de placas impresas; extender comparación a Tet10 y láminas. | Reproducir cada probeta y su montaje como caso separado, conservar adimensionalización y auditoría de réplicas. No atribuir sus propiedades al CAD. |
| [Taguchi et al. 2020](https://www.jstage.jst.go.jp/article/tjsass/63/1/63_T-19-27/_article) | Referencia de morphing pasivo, forma deformada, CL y Cp; zona flexible posterior. | Re≈1,78×10⁵, relación geométrica y condiciones bidimensionales del ensayo deben distinguirse del álabe finito. |
| [Kai et al. 2023](https://doi.org/10.1007/978-981-19-2689-1_78) | Arquitectura impresa con mecanismo compliant; orientar representación de piel y estructura. | Extraer constitutiva, geometría, uniones y montaje antes de crear una réplica cuantitativa. |
| [Ochi et al. 2024](https://doi.org/10.2514/1.J063000) | Referencia de acoplamiento estático y efecto de elegir aerodinámica viscosa. | El resumen disponible no proporciona todas las entradas; obtenerlas antes de usar curvas como objetivo numérico. |
| [Malim et al. 2023](https://doi.org/10.3390/aerospace10030256) | Propiedades AM iniciales Tough PLA/ABS que ya utiliza el código. | Registrar impresión/orientación; Poisson 0,35 del ejemplo es supuesto y no identifica amortiguamiento. |
| [Bathe y Dvorkin 1985](https://doi.org/10.1002/nme.1620210213) | Base de verificación de placa y decisión informada al desarrollar láminas. | Verificar por separado la formulación curva y sus conexiones al sólido. |
| [SHARPy, del Carre et al. 2019](https://joss.theoj.org/papers/10.21105/joss.01885) | Referencia independiente de arquitectura FSI y UVLM, incluyendo estela y dinámica. | Igualar hipótesis y operadores en benchmarks; su estructura de viga no reemplaza directamente el CAD interno. |
| [Ōtomo et al. 2025](https://researchportal.bath.ac.uk/en/publications/scaling-unsteady-load-alleviation-in-airfoils-with-flexible-trail/) | Referencia posterior para pieles flexibles y escalas de deformación/carga no estacionaria. | El ensayo tiene movimiento de plunging y condiciones hidrodinámicas: no es validación directa de un álabe fijo en aire. |

Para cada trabajo guardar una ficha con DOI/versión, geometría, material,
montaje, fluido, Reynolds, excitación, señales disponibles, normalización,
incertidumbre, licencia y páginas/figuras fuente. Si solo hay curvas, registrar
digitalización y su error. Separar mediciones de otros resultados numéricos.
Mantener tres clases: benchmark reproducido, comparación de tendencia y
referencia metodológica. No agruparlas como si aportaran la misma validación.

## Organización propuesta del código y los resultados

Reutilizar las interfaces existentes y añadir capacidades mediante PR pequeños:

| Pieza propuesta | Responsabilidad |
| --- | --- |
| `cases/cad_original` o contrato equivalente | Configuración física y numérica versionada, sin rutas personales obligatorias. |
| Preprocesamiento de regiones CAD | Etiquetas físicas, espesores, superficies mojadas, sensores y calidad de malla. |
| `fem3d.tet10` y soporte asociado | Elemento, cuadratura, frontera Tri6, ensamblaje, masa, tensiones y cargas consistentes. |
| Operador de presión sobre piel curva | Distribución física de cargas y transferencia adjunta verificada. |
| Restricciones de mecanismo | Guías, uniones y posteriormente contacto/fricción documentados. |
| Observación del CAD | H/O y sensores de cola comunes a mallas y experimentos. |
| Ejecutor de campaña CAD | Casos estáticos, modales, forzados, refinaciones y comparaciones. |
| Aerodinámica no estacionaria | Estado de circulación/estela separado del estado estructural. |

Estos nombres son propuestas, no funciones disponibles. Hoy se reutilizan
`fem3d.read_mesh`, `assemble`, `traction`, `modal`, `newmark`, `transfer`,
`one_way`, `pae.vlm` y las funciones `ael` documentadas en la API.

Cada corrida guardará configuración, hash/commit, regiones, malla/calidad,
versiones, carga aplicada, balances, convergencias, sensores, señales y alcance.
Las salidas permanecerán en `results/`, excluido de Git. Registrar si geometría,
material, presión o amortiguamiento son supuestos. Datos externos conservan
su licencia y procedencia. No publicar el CAD original por defecto.

Las figuras principales serán: forma deformada CAD con escala indicada,
w_TE frente a envergadura, modos relevantes, señal temporal y espectro/FRF,
comparación experimental con incertidumbre y, cuando exista solver adecuado,
estela con magnitud física identificada. Evitar figuras decorativas de flujo.

## Matriz de casos y presupuesto de cálculo

Construir campañas pequeñas antes de ampliar los barridos:

| Grupo | Casos mínimos | Qué diferencia permite interpretar |
| --- | --- | --- |
| Geometría y discretización | CAD original, tres mallas y carga analítica idéntica. | Error numérico sin mezclar cambios de carga o material. |
| Mecánica | Sin flujo, varios niveles de carga/descarga y configuraciones identificadas del slider. | Flexibilidad, no linealidad, histéresis y efecto del montaje. |
| Aerodinámica rígida | Forma CAD inicial inmóvil, varias incidencias y Reynolds de referencia. | Error de cargas antes de atribuirlo a aeroelasticidad. |
| Deformación desacoplada | Misma forma inicial y presión calculada/medida, estructura flexible. | Respuesta estructural para un campo de presión conocido. |
| Efecto de forma | Fluido sobre la forma deformada congelada, sin movimiento adicional. | Cambio aerodinámico por geometría frente a realimentación dinámica. |
| Acoplamiento | Una vía y dos vías, iguales condiciones iniciales y físicas. | Magnitud de la realimentación y error del procedimiento desacoplado. |
| Dinámica | Perturbación inicial, presión prescrita conocida y después flujo acoplado. | Separar vibración libre, forzada e interacción fluido–estructura. |

Elegir Uref y los ángulos con el equipo experimental. Una continuación
computacional inicial puede usar fracciones 0,25/0,5/0,75/1 de Uref, con
incidencias pequeñas dentro del régimen comprobado y refinación donde cambie
la respuesta. No extrapolar ese régimen hasta pérdida con VLM. El valor
histórico de 10 m/s no constituye un límite ni una condición experimental confirmada.

Para cada configuración, variar primero una entrada a la vez: espesor, módulo,
rigidez del apoyo, posición del slider y amortiguamiento. Después construir
un conjunto combinado de incertidumbre si las mediciones lo justifican.
Conservar controles rígidos y las mismas escalas en todas las comparaciones.
No asumir que mayor flexibilidad debe aumentar siempre CL o reducir vibración.

Antes de resolver todo el CAD con Tet10, contar aristas únicas y grados de
libertad, estimar matrices dispersas y medir ensamblaje/factorización en un
submodelo representativo. Tet10 tiene matrices locales 30×30 frente a 12×12
de Tet4: el almacenamiento local y el llenado de la factorización requieren
un presupuesto específico, no una extrapolación por número de elementos.
Ensamblar por bloques y reutilizar operadores/factores cuando sea válido.
Registrar memoria máxima y tiempo junto a la precisión obtenida.

Si el coste supera los recursos disponibles, priorizar remallado local,
submodelos de flexuras con condiciones de interfaz verificadas, reducción
modal y eventualmente un híbrido sólido–lámina derivado del mismo CAD. La
decisión debe conservar masa, conectividad y flexibilidad demostradas. Definir
recursos de CPU/RAM antes de comprometer campañas 3D viscosas o grandes barridos;
no asignar fechas de entrega sin medir primero el coste del caso piloto.

## Primer bloque de trabajo ejecutable

Este es el orden recomendado para retomar implementación sin empezar por FSI:

1. Congelar la corrida histórica y crear la configuración física del CAD.
2. Localizar pieles, borde de salida, flexuras, zonas esbeltas y apoyo real.
3. Definir 5–9 estaciones TE y medir desplazamiento firmado, no solo máximo global.
4. Aplicar una presión analítica en un parche fijo sobre piel identificada.
5. Comparar ambas mallas con esa carga, apoyo y sensores idénticos.
6. Separar el efecto de huella de carga del efecto de elemento y geometría.
7. Implementar Tet10 y Tri6, con benchmarks de flexión, masa y modos.
8. Generar tres mallas locales del CAD y comparar w_TE, energía y frecuencias.
9. Identificar apoyo/material/slider con ensayos o mantener escenarios explícitos.
10. Solo entonces aplicar presión aerodinámica inicial y temporal al CAD convergido.

El primer hito útil es una curva presión–deformación y una familia modal del
CAD con error numérico acotado. Este hito es valioso aunque aún falten
presiones experimentales o un solver de estela.

## Cómo retomar el proyecto

Repositorio: [passive-airfoil-aeroelasticity-matlab](https://github.com/MC-Andes/passive-airfoil-aeroelasticity-matlab).
El STEP local del usuario está en
`/Users/nuseche/Downloads/Airfoil with slider alone.stp`; en otra máquina se
debe configurar su ubicación y comprobar el hash antes de regenerar mallas.
Las mallas publicadas no requieren ese archivo para ejecutarse.

Desde la raíz, con MATLAB R2026a:

```matlab
addpath('matlab')
run_verification();
run('matlab/inspect_reference.m');
run('matlab/run_fem_reference.m');
run('matlab/run_fem_convergence.m');
run_tail_benchmarks();
```

La última orden compara placas de literatura; no ejecuta un modelo de placa
del mecanismo original. Los scripts de referencia sobrescriben resultados
homónimos: guardar una copia de la corrida histórica antes de repetirlos.
La malla fina consume más tiempo/memoria. Utilizar Cholesky disperso con
ordenamiento, como ya hace el solver; evitar factorizaciones densas.

En terminal, `python3 tools/check_repository.py` verifica integridad. Para
un documento nuevo todavía sin seguimiento Git, añadirlo al índice o comprobar
explícitamente sus enlaces, porque la herramienta inspecciona archivos versionados.
No ejecutar herramientas Python de CAD sin preparar sus dependencias opcionales;
el núcleo numérico seguirá desarrollándose y verificándose en MATLAB.

Al terminar cada bloque actualizar: estado de la etapa, archivos y commit,
pruebas/comparaciones ejecutadas, supuestos pendientes, criterio de cierre y
siguiente acción concreta. Mantener este plan como referencia de aceptación,
sin confundir tareas planificadas con capacidades implementadas.


## Ejecución y criterios de cierre — continuación del 17 de octubre de 2026

El primer bloque estructural fue rechazado: se conserva su
[registro](plan-cad-first-block.md) y sus [resultados](cad-original-history.md).
Se continuó el trabajo sobre la misma geometría, sin cambiar el objetivo ni
ajustar E para compensar mallas deficientes.

El usuario confirmó que no había datos de montaje, material, slider o ensayos,
y autorizó: «Ninguna, suponla a traves del analisis de la pieza real y ejemplos
de literatura». Esa autorización permite escenarios **condicionales**, con
hipótesis trazables; no crea mediciones ni valida físicamente el espécimen.
El [contrato de hipótesis](../cases/cad_original/assumptions.json) y las
[fichas de literatura](../data/benchmarks/cad_literature.json) fijan esa distinción.

| Etapa | Desarrollo ejecutado | Evidencia / alcance de aceptación |
| --- | --- | --- |
| 1 Referencia | Respaldo histórico con hashes, reproducción, contratos físicos/numerales y fuentes exactas por campaña. | Diferencia histórica de desplazamiento menor de 1,3e−8; siete suites MATLAB, 69 controles pasando. |
| 2 Geometría | Cap original con todas sus curvas, 587 firmas, espesores normales, ocho puertos de cavidades, barra candidata 7/9, nueve sensores y reflexión X/Z/Y. | Grupos y escenarios ligados al SHA del STEP; presión interna ambiente como hipótesis de ventilación, con análisis de sensibilidad de impedancia. |
| 3 Transferencia | Parches físicos Tri3/Tri6, inversión curva, H/H', fuerzas de VLM con su brazo correcto y tracciones viscosas sobre ambas pieles. | Fuerza/momento/trabajo verificados; benchmark analítico de paraboloide y comparación de huella histórica. |
| 4 Elementos | Tet10 isoparamétrico, integración por bloques, masa, rigidez, frontera, apoyos, tensiones y geometría proyectada al CAD. | Tres niveles: último cambio TE 0,270%, energía 0,526%, f1 0,154%; comparación FEM curva completa: TE 0,00190%, f1 0,00620%. |
| 5 Mecánica | A monolítico, B bloqueado, C guiado; multiplicadores, contacto/fricción incremental, gravedad, sólido total Lagrangiano y sensibilidades. | A/B/C a 100 Pa: 0,22424/0,01234/0,06794 mm; tres mallas. Guías/material/apoyos son hipótesis, no identificación experimental. |
| 6 Modos/dinámica | Modos completos del CAD, MAC sobre cuerda/envergadura, energía modal regional, FRF, Newmark completo/reducido y corrección estática residual. | Banda A 0–60 Hz; f1=20,591 Hz. Guías B/C f1=58,222/38,531 Hz. Error de sensor FEM completo/reducido 5,46e−7 en el primer subespacio; truncamiento evaluado por separado. |
| 7 Una vía | Superficie media desde pieles reales, refinación de paneles/herraduras, fuerzas físicas conservativas y solución FEM completa con carga refinada. | Con 512×26, error FEM completo/reducido ≤0,301% en cargas rígidas y equilibradas; tracción viscosa de sección 0,836%. Se guardan presión/movimiento sincronizados. Cp de túnel no se atribuye al VLM. |
| 8 Dos vías estático | Actualización de geometría/normales, continuación, residuos, límites de rigidez/carga y sensibilidad de relajación. | Mapas A/B/C y dos idealizaciones estacionarias separadas. Último cambio TE: herraduras 0,413%, anillos 0,539%; discrepancia entre modelos declarada. Límites de carga/rigidez y relajación comprobados. |
| 9 UVLM | Anillos, Kutta/Kelvin, Bernoulli no estacionario, estela congelada/libre, convección y derivadas/transferencia verificadas. | Theodorsen/Sears en k=0,5/1,6/6,4: menos de 2% amplitud y 2° fase; movimiento prescrito CAD: error de fuerza 5,81e−5. Subdivisión temporal exacta contrastada con ensamblaje directo, sin comprimir historia. |
| 10 FSI dinámico | Acoplamiento fuerte, eliminación exacta de circulación, historial completo, perturbación libre, modos/dt/paneles, velocidades y bandas materiales. | Intercambio no lineal CAD/tangente 0,0241%; refinación 96/128: decaimiento 1,63%. Paso final A 31,25→15,625 µs: 1,08%, RMS 0,54%, fase 0,41°. Tiempo/espacio por velocidad y guías B/C aceptados; bandas condicionales calculadas. |
| 11 Viscoso | Navier–Stokes propio 2D/3D, presión/Brinkman acoplados, masa/cantidad de movimiento, pared móvil, pieles y estela calculada. | Poiseuille/Taylor–Green/ABC 3D y cilindro Re20 Cd≈2,0475 frente a 2,05. Casos CAD Re20/100, perfiles de estela con cambio <0,5%; Cl conserva dispersión 4–5% y tracción de piel error adicional. Presupuesto de túnel: unos 10 TB de trabajo. No se afirma DNS/URANS/LES validado a Re de túnel. |

El [informe principal](cad-original-results.md) concentra números, carpetas,
figuras, refinaciones y límites. La [guía estructural](cad-original.md), la
[ruta acoplada](cad-coupled.md) y la [rama viscosa](viscous-flow.md) permiten
repetir los componentes. El código y las entradas públicas siguen separados
del STEP y de los resultados locales, que permanecen fuera de Git.

El cierre algebraico de la nueva campaña usa el estimador dual de energía
que se declaró antes de ejecutarla; conserva los residuos directos y explica
sus límites en [residuos numéricos](numerical-residuals.md). La sexta frecuencia
no cumple 1% entre las últimas mallas y no se presenta como convergida. Una
norma estimada y una diferencia pequeña entre mallas no son cotas certificadas.

La validación experimental del material/montaje, Cp de túnel, flutter real y
PIV no se han fabricado ni se sustituyen por las hipótesis autorizadas. El
resultado del desarrollo se interpreta como estudio numérico reproducible y
condicional del álabe original, con las extensiones y límites físicos explícitos.


### Cierre registrado de la ejecución numérica

La revisión conjunta acepta el modelo condicional: `acceptedConditionalNumericalModel=true`
en `results/cad_original/final_review/acceptance.json`. El manifiesto de entradas,
las tablas y los hashes permiten revisar cada criterio sin seleccionar carpetas
por fecha. El informe local con figuras está en
`results/cad_original/final_review/report.md`; el informe metodológico público
se mantiene en [resultados del CAD](cad-original-results.md).

Las once etapas tienen desarrollo y evidencia numérica dentro de sus hipótesis.
La comprobación Duffy 4/6 de todos los elementos curvos, en sus seis modos
y campo estático, cambia las energías menos de 8,27e−10 relativo; conserva
por separado los controles de jacobianos y la comparación geométrica completa.
Se calcularon también forma posterior, pendiente y curvatura en líneas materiales:
la curvatura ajustada cambia 3,20% entre las dos últimas mallas y se informa
como observable secundario, sin extenderle la aceptación del 2% de TE/energía.
Tampoco se atribuye esa precisión a Cl viscoso, tensiones singulares o Cp de piel.

Las respuestas libres A/B/C decaen bajo las hipótesis estudiadas; no se necesitó
introducir una presión oscilante para producir el transitorio. No se establece
un umbral experimental de flutter ni una oscilación autosostenida real. Las
mediciones de material/montaje/slider, Cp y PIV siguen siendo evidencia externa
no disponible; la autorización de hipótesis permitió completar el estudio
condicional, sin fabricar esa validación.
