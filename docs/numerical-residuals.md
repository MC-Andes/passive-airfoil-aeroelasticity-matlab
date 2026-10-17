# Residuos y precisión en sólidos esbeltos

El CAD contiene aristas submicrométricas y tetraedros con razones geométricas
cercanas a 10⁻⁶. La resta de términos grandes en `K*u` puede dominar un
residuo normalizado por una fuerza aerodinámica pequeña. Una factorización
estable no elimina la sensibilidad de ese residuo a la representación.

La evaluación adicional usa exactamente la misma elasticidad y cuadratura,
aplicando gradientes, deformación y tensión por elemento. Se eliminan los seis
movimientos rígidos locales por una proyección ortogonal antes de diferenciar
y después de integrar fuerzas. En aritmética exacta, B anula esos movimientos
y el operador no cambia; la proyección controla cancelación numérica. No añade
resortes, masa ni estabilización mecánica. Se contrasta con la matriz ensamblada
sobre vectores aleatorios, deformaciones afines y movimientos rígidos.

Se mantienen los residuos de fuerza y del problema propio en sus unidades y
normalizaciones anteriores. También se calcula el residuo dual de energía:

`eta_E = sqrt(r' * Kff^{-1} * r / (u' * K * u))`.

Para el problema lineal SPD, el numerador es la energía del error de solución
si el inverso se aplica exactamente. Aquí se estima con el factor disperso y
se etiqueta como estimación a posteriori, sin atribuirle una certificación de
intervalo. El residuo modal usa el mismo producto dual, normalizado por la
energía modal. Este producto mide la proximidad espectral de 1/lambda en la
transformación auto-adjunta por K, y evita interpretar ruido de alta rigidez
como error del desplazamiento flexible dominante.

Antes de ejecutar la campaña fina se adoptó 10⁻⁶ para esta medida energética,
conservando todos los residuos de fuerza originales y los balances globales.
Además se mantienen los criterios físicos de refinación del 2% en sensores y
energía, y del 1% en frecuencias de la banda declarada. No se ajusta E para
compensar discretización ni se omiten las divergencias entre mallas.

La banda inicial asumida es 0–60 Hz: contiene el primer modo de cola, cercano
a 20 Hz, y la escala convectiva U/c≈40 Hz a U=10 m/s. Se reportan también modos
superiores y su error, pero no se atribuye resolución validada fuera de esa
banda sin ampliarla y refinar. Esta elección es una hipótesis de estudio ante
la ausencia de datos, autorizada por el usuario; no es ancho de banda medido.
