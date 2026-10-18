# Registro histórico de resultados del bloque CAD original — 17 de octubre de 2026

El bloque de software se implementó y ejecutó sobre el CAD original. **El modelo estructural todavía no cumple los criterios de cierre**: no se declara una predicción validada de la cola, del slider o de su vibración en aire. El [plan](plan-cad-original.md) conserva sus objetivos y dependencias. La [guía ejecutable](cad-original.md) describe contratos y comandos.

## Configuración y procedencia

- STEP sin modificar, SHA256 `5bacb470df1b1239a0d8af2181f10ca14077a4c0eba3b3531126c8981ea4793e`.
- Escenario A: sólido monolítico, Tough PLA isotrópico publicado como referencia y empotramiento supuesto en x/c≤0,30. No se identificó el montaje real ni el slider.
- Presión analítica impuesta de 0/25/50/100 Pa sobre piel superior expuesta, x/c=0,78..0,84 y eta=0,25..0,75. No procede de VLM ni de una medición.
- Nueve sensores materiales en x/c=0,999, eta=0,05..0,95. Se informa desplazamiento z firmado; las magnitudes de la tabla son máximos absolutos de esos sensores.
- MATLAB R2026a Update 4, host con 24 GiB de RAM. Base Git `9d77047876ba6705589eac924fccf1424606448c`, cambios locales sin commit; cada campaña registra hashes de fuentes e inputs.
- Resultados nuevos bajo `results/cad_original/`; los datos publicados y el STEP permanecen sin modificaciones.

## Campañas ejecutadas

| Elemento y nivel | Nodos | Elementos | Máximo de sensores a 100 Pa [mm] | Residuo estático relativo |
| --- | ---: | ---: | ---: | ---: |
| tet4_level_1 | 21,149 | 69,082 | 0.006653 | 6.69e-08 |
| tet4_level_2 | 32,925 | 109,186 | 0.018932 | 3.14e-07 |
| tet4_level_3 | 85,367 | 286,037 | 0.018284 | 4.98e-07 |
| tet10_level_1 | 112,829 | 58,097 | 0.176759 | 1.43e-05 |
| tet10_level_2 | 131,612 | 69,082 | 0.169009 | 2.8e-05 |
| tet10_level_3 | 206,099 | 109,186 | 0.183355 | 6.38e-05 |

Las tres mallas Tet4 tienen 69 082, 109 186 y 286 037 elementos. Las tres Tet10 corresponden a remallados nominales de 20, 15 y 12 mm, con aristas rectas. Las dos familias no constituyen el mismo conjunto de tamaños: la comparación directa de orden utiliza el Tet4 de 69 082 elementos y el Tet10 del mismo mallado, con geometría y caras idénticas.

| Criterio entre los últimos niveles | Tet4 | Tet10 | Objetivo del plan |
| --- | ---: | ---: | ---: |
| Cambio del vector de sensores | 5.375% | 7.824% | <2% |
| Cambio de energía | 0.543% | 8.678% | <2% |
| Cambio de primera frecuencia emparejada | 2.842% | 3.643% | <1% |

El cambio de sensores es la norma infinito de la diferencia, dividida por la norma infinito del nivel final; no equivale al cambio de un único máximo. La secuencia no es monótona y varios modos tienen emparejamiento insuficiente en esos sensores. El MAC incluye xyz y pesos de envergadura; los modos próximos o poco observables requieren revisar subespacios/sensores y no se consideran convergidos por su número.

Los residuos máximos del problema modal Tet10 llegan aproximadamente a 1,34×10⁻⁴. El equilibrado y refinamiento no eliminan el límite numérico de estas mallas esbeltas. Las soluciones, curvas de presión y FRF se guardan como **diagnósticos condicionales**. No hay una cota aceptada de error para extrapolar una predicción física.

La elevación Tet4→Tet10 sobre la misma malla gruesa cambia el máximo de sensores de 0,006653 a 0,169009 mm a 100 Pa. Los benchmarks independientes pasan, pero esta diferencia no demuestra que el segundo valor sea exacto. Primero hay que mejorar y volver a verificar las mallas locales y la representación de detalles del CAD.

## Transferencia y geometría

El área proyectada del parche es 0,002667558534 m² en todos los niveles: a 100 Pa, la componente vertical integrada es −0,2667558534 N. Las áreas reales son de aproximadamente 0,00270154 m² y se registran por nivel. Los defectos de fuerza, momento y trabajo cumplen 10⁻⁶ en las pruebas controladas; la presión no se aplica a cavidades por proximidad.

Con las mismas fuerzas de integración y nodos de la misma piel, el operador histórico de 16 vecinos cambia el vector de sensores un 0,0430% en la gruesa y un 0,0133% en la fina. Su radio máximo pasa de 21,70 a 9,70 mm. Esta comparación aislada no explica la diferencia estructural observada. No valida la carga VLM histórica sobre superficies medias ni su brazo de aplicación.

Cada malla tiene atlas con 587 firmas de cara, coordenadas de los 200 peores tetraedros e intervalos verticales de material. En la malla gruesa, el elemento peor tiene razón singular 1,55×10⁻⁶ cerca de (x,y,z)=(0,219748; 0,018752; 0,057417) m canónicos. El atlas permite localizar el defecto; no se le asigna una función de flexura o slider sin evidencia. Los intervalos verticales no son espesor normal de piel.

Los intentos de 10 mm y de optimización Netgen no se usaron como resultados: uno produjo fallo nativo y otro malla no válida. Se conservaron logs y se rechazaron conectividades degeneradas/no manifold. Las mallas válidas adicionales se obtuvieron con 12 y 20 mm, sin borrar caras pequeñas o uniones.

## Recursos y reproducción histórica

La campaña Tet10 en serie tardó unos 140 s en este host, incluyendo arranque y exportaciones: `/usr/bin/time -l` registró RSS máximo de 10,88 GB y `peak memory footprint` de 15,41 GB. Son métricas distintas del sistema operativo; no equivalen a los bytes de matrices. El nivel de 618 297 DOFs se midió, mientras que elevar la malla publicada fina produciría 1 602 759 DOFs y sigue fuera del presupuesto actual. No se promete ese coste por extrapolación lineal.

La corrida histórica se respaldó antes de repetirla, con hashes en `history_9d77047/snapshot_manifest.json`. Sus scripts se ejecutaron después con idénticas entradas y salida redirigida a `historical_reproduction/`; el equilibrado nuevo reproduce 0,00552013494 y 0,01628631883 mm de desplazamiento máximo. La diferencia relativa máxima respecto al respaldo es menor de 1,3×10⁻⁸. El cambio histórico de 66,1% entre mallas continúa sin converger. Además, regenerar la malla de 20 mm con las opciones explícitas y sin preferencias personales produce CSV de nodos, elementos y caras idénticos byte por byte; la comprobación se guarda en `review_20261017/reproducibility.json`.

## Pruebas y alcance

`run_verification()` pasa 44 controles: 33 históricos y 11 nuevos. Tet10 verifica masa, rigidez, solución afín, permutaciones, flexión y frecuencia analíticas, cuadratura curva, geometría plegada, Tri6, observación, transferencia, vibración libre y procedencia. El chequeo de integridad de inputs y documentación también pasa con `--include-working-tree`.

Las FRF usan amortiguamiento modal supuesto del 1%. El decaimiento libre adicional impone 10 μm de perturbación inicial, sin presión de flujo; compara dos pasos temporales con la solución analítica amortiguada y verifica energía. El error de amplitud baja de 1,092% a 0,273%, y el defecto energético de 0,0706% a 0,0176%. Su convergencia temporal no valida los modos espaciales ni el amortiguamiento físico.

## Estado y siguiente trabajo

1. Mantener pendientes los criterios estructurales: localizar/remallar detalles esbeltos y encuentros sin cambiar conectividad física; estudiar tamaño por espesor/curvatura, calidad del jacobiano y geometría curva Tet10. Repetir al menos tres niveles con estas mismas cargas y sensores.
2. Identificar montaje, material/proceso, masa, espesores y mecanismo del slider con el equipo experimental. Los escenarios B/C siguen bloqueados por falta de interfaces y cinemática documentadas.
3. Resolver observación modal de modos cercanos/locales, tolerancias algebraicas y error de reducción en la banda objetivo. No ajustar E para compensar la falta de convergencia.
4. Tras cerrar estructura y cargas, continuar con aerodinámica fija, respuesta desacoplada y las etapas 8–10. UVLM y estela viscosa no se han implementado en este bloque ni se presentan como capacidades existentes.

Campañas aceptadas como ejecución:

- Tet4: `results/cad_original/campaign_20261017T182550399`.
- Tet10: `results/cad_original/campaign_20261017T182259762`.
- Referencia histórica reproducida: `results/cad_original/historical_reproduction/`.
- Comparación gráfica y procedencia de posprocesado: `results/cad_original/review_20261017/`.

Los JSON de aceptación mantienen `structuralModelClosed=false`; ambos cierres de discretización son falsos y los cierres algebraicos también son falsos (en Tet4 fallan residuos modales, y en Tet10 también los estáticos). Este estado impide declarar terminadas las once etapas.


Este registro describe el primer bloque rechazado. El estado actual está en el
[informe principal](cad-original-results.md).
