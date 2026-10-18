# Registro del primer bloque estructural

Este bloque histórico fue rechazado y queda superado por el
[estado actual del plan](plan-cad-original.md).

## Ejecución registrada — 17 de octubre de 2026

Se implementó y ejecutó el bloque estructural del CAD, manteniendo el objetivo
original y el núcleo MATLAB propio. El [informe de resultados](cad-original-results.md)
contiene cifras y limitaciones; la [guía de campaña](cad-original.md) permite
repetirlo. El trabajo está en el árbol local basado en `9d77047`, sin crear un
commit o PR nuevo. Los resultados y el CAD original siguen fuera de Git.

**Estado de aceptación: no cerrado.** Se ejecutaron tres mallas Tet4 y tres
Tet10, pero ninguna familia cumple los criterios estructurales. Entre los
últimos niveles cambian los sensores 5,375%/7,824% y la energía 0,543%/8,678%
(Tet4/Tet10). La primera frecuencia emparejada cambia 2,841%/3,643%.
Los residuos modales también exceden el criterio algebraico en parte de las
mallas. No se modificaron tolerancias para declarar éxito.

| Etapa | Estado comprobado tras este bloque | Criterio pendiente / siguiente acción |
| --- | --- | --- |
| 1 Referencia | Cerrada para la reproducción numérica disponible: respaldo con hashes, contrato físico/numerical, fuentes e inputs trazados y 44 verificaciones. La referencia histórica se reproduce con diferencia relativa de desplazamiento menor de 1,3×10⁻⁸. | Material y montaje siguen clasificados como supuestos; la reproducción no los valida. |
| 2 Geometría | Parcial: 587 firmas por malla, pieles candidatas, exposición del parche, intervalos verticales de material, nueve sensores y transformación axial verificada. | Identificar función mecánica de nervaduras/flexuras/slider, apoyo real, espesores normales y superficies curvas con evidencia. |
| 3 Transferencia | Verificada para el parche analítico fijo: recorte físico, conservación y comparación aislada con 16 vecinos. Integración Tri6 curva verificada por separado. | Completar la transferencia de cargas aerodinámicas en sus superficies reales; no dar por resuelto el brazo de carga del VLM. El recorte físico curvo exige integración paramétrica adicional. |
| 4 Elementos | Tet10/Tri6 implementados y verificados, con rigidez, masa, frontera, apoyos, conectividad, tensiones de centroide y posprocesado. Tres niveles de cada familia ejecutados. | Cierre rechazado: remallar detalles esbeltos y encuentros, controlar espesor/curvatura, jacobianos y error algebraico; repetir tres niveles. No se proyectaron nodos al CAD ni se borraron uniones. |
| 5 Mecánica | Escenario A explícito; B/C se rechazan si se solicitan sin interfaces/cinemática identificadas. | Ensayos de material, apoyo y slider; restricciones de mecanismo, contacto y no linealidad permanecen pendientes. |
| 6 Modos/dinámica | Seis modos por nivel, MAC en sensores comunes, error estático por truncamiento, FRF con amortiguamiento supuesto y decaimiento libre con refinación temporal/energía. | Convergencia espacial y residuos modales; modos locales/próximos y subespacios; amortiguamiento medido y comparación completa/reducida en la banda objetivo. |
| 7 Una vía CAD | Curvas 0/25/50/100 Pa y forma deformada bajo presión analítica sobre la piel exterior. | Cerrar estructura y aerodinámica fija del CAD; usar después presiones calculadas o medidas con procedencia. Las cargas analíticas no son presión de túnel. |
| 8 Estático dos vías | Pendiente, sin anticipar cierre estructural. | Cumplir dependencias 3–7 antes de implementar/evaluar equilibrios como predicciones. |
| 9 UVLM | Pendiente; el VLM estacionario continúa siendo la única aerodinámica de esta ruta. | Implementación propia y benchmarks independientes de circulación, fuerzas no estacionarias y estela. |
| 10 FSI dinámico | Pendiente. | Fluido y estructura verificados/convergidos, amortiguamiento identificado y balances de intercambio. |
| 11 Estela viscosa | Pendiente, como extensión física específica. | Definir pregunta y datos de estela, presupuesto y solver viscoso verificado. No se reemplazó UVLM por líneas decorativas ni por un motor externo. |

### Archivos y pruebas del bloque

- Contratos: `cases/cad_original/case.json` y `tet10-convergence.json`.
- Ejecución: `run_cad_original`, `run_cad_modal_response` y paquete `cad`.
- Sólido: `fem3d.tet10`, funciones de forma, integración superficial,
  elevación conforme y adaptación de los operadores FEM existentes.
- Solución: Cholesky disperso equilibrado, refinamiento de residuos con el
  operador original y modos calculados con congruencia de K/M.
- Pruebas: los 33 controles históricos más 11 de CAD/Tet10, todos pasando
  en MATLAB R2026a Update 4; integridad de los 19 inputs publicados y enlaces
  comprobada con `python3 tools/check_repository.py --include-working-tree`.
- Preprocesamiento: rechazo de elementos degenerados y fronteras no válidas,
  sin reparación oculta de topología; un intento de 10 mm y Netgen fallidos
  quedaron excluidos y registrados.
- Recursos: el mayor Tet10 ejecutado tiene 618 297 DOFs. La campaña en serie
  registró unos 140 s, RSS máximo de 10,88 GB y footprint máximo de 15,41 GB
  en un host de 24 GiB. No se extrapoló ese coste a 1,6 millones de DOFs.

Salidas principales: `results/cad_original/campaign_20261017T182550399`
(Tet4), `campaign_20261017T182259762` (Tet10),
`historical_reproduction` (referencia reproducida) y `review_20261017`
(comparación gráfica y hashes de posprocesado), todas bajo
`results/cad_original/`. El decaimiento modal condicional está en
`campaign_20261017T182259762/tet10_level_3/free_decay`.

### Siguiente acción concreta y supuestos pendientes

Usar las coordenadas de `atlas/worst_elements.csv` y los intervalos de material
para preparar un remallado local de flexuras/encuentros, con controles de
espesor, curvatura y calidad. Mantener exactamente el parche, sensores y
apoyo del contrato al compararlo. La falta de convergencia no se corrige
ajustando E ni cambiando el mecanismo para conseguir una respuesta deseada.
Si el sólido completo supera los recursos, verificar submodelos derivados
del mismo CAD o la alternativa híbrida prevista en la etapa 4.

En paralelo, obtener del equipo experimental la definición del montaje,
material/proceso, masa/espesores, operación del slider y amortiguamiento.
Sin esa evidencia se conserva A como hipótesis y no se habilitan B/C.
El informe registra una primera familia de curvas/modos útil para diagnóstico;
el hito de predicción con error numérico acotado y las once etapas completas
**siguen pendientes**, conforme a sus criterios originales.
