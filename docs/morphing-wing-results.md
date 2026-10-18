# Comprobaciones ejecutadas con el CAD actualizado

Esta evidencia corresponde exclusivamente al STEP
`a5b289434de68e737d9c0bed552d32c642f6be22fc84d8c71d49434c501f14da`, PETG
isótropo E=1,6 GPa, ν=0,35 supuesto y ρ=1.270 kg/m³, apoyo en el alojamiento
adherido y frente del slider rígido idealizado. La inserción 0,324 m es una
recomendación asumida, no un dato final. Véase la
[definición completa](morphing-wing.md).

## Geometría y malla nuevas

La auditoría adaptativa encuentra un sólido BREP válido, 630 caras y volumen
0,0007090739570 m³. El incremento de volumen frente al CAD anterior es
15,1976%. La geometría no es prismática. Las caras 71–74 definen el alojamiento
cuadrado de 15,2 mm; las caras 7–9 forman el frente de accionamiento. Las
firmas se contrastan con el STEP antes de usar etiquetas. La diferencia de
0,2 mm entre alojamiento y eje nominal se interpreta como holgura de adhesivo;
la capa no se representa como material independiente.

Se generó una malla 3D de 68.453 nodos Tet4, 224.188 tetraedros y 129.766 caras
de frontera. Al elevar a Tet10 conserva los tetraedros y alcanza 425.987 nodos,
1.277.961 grados de libertad. El error de volumen es 0,1001%. La topología y
la correspondencia frontera/volumen pasan; hay 14.753 elementos con razón de
valores singulares del jacobiano menor de 1e−4. El determinante normalizado
mínimo es 1,56e−9: son detalles muy esbeltos, no una malla demostrada óptima.
No se borran elementos o caras para mejorar las cifras.

Fuentes: `results/morphing_wing/geometry/`,
`results/morphing_wing/meshes/pilot/` y
`data/reference/morphing_wing_geometry.json`.

## Prueba tangente, desplazamiento del slider −1 µm

Se resuelve el sólido completo, sin presión ni fuerza impuesta, con el comando
como condición esencial del slider. Las reacciones equilibran el apoyo y el
trabajo de frontera coincide con dos veces la energía elástica lineal. Es una
comprobación pequeña alrededor de slider cero, no un cálculo a recorrido final.

| Magnitud | Tet4 | Tet10 |
| --- | ---: | ---: |
| DOFs | 205.359 | 1.277.961 |
| Masa [kg] | 0,899622 | 0,899622 |
| Máximo absoluto TE [µm] | 4,19703 | 4,26746 |
| Reacción del accionamiento en x [N] | −0,164510 | −0,021166 |
| Energía [J] | 8,2255e−8 | 1,0583e−8 |
| Frecuencia fundamental, slider bloqueado [Hz] | 180,923 | 45,1251 |
| Residuo estimado de energía estática | 7,97e−9 | 3,12e−7 |
| Máximo von Mises de centroide [MPa] | 0,50006 | 0,021167 |

El vector TE cambia aproximadamente 1,65% al cambiar de elemento, pero la
reacción/energía de Tet4 es **7,77 veces** la de Tet10 y su primera frecuencia
**4,01 veces** mayor. Por tanto, el acuerdo del desplazamiento de salida no
justifica aceptar la rigidez, las tensiones o los modos de Tet4. Se elige Tet10
como referencia y **no se declara convergencia mecánica del modelo nuevo**:
faltan refinamientos con el mismo elemento, condiciones físicas y observables.
La comparación actual es de orden de elemento sobre una sola malla base.

Las seis frecuencias Tet10 del ensayo condicional con slider bloqueado son
45,1251; 93,9173; 129,7505; 141,1559; 142,0105; 142,4492 Hz. No equivalen a
frecuencias medidas. Los residuos modales estimados de energía quedan por
debajo de 4e−8. El máximo de tensión es local, dependiente de restricciones,
geometría y malla; no es un criterio de rotura FDM ni un factor de seguridad.

Corridas preservadas:

- Tet4: `results/morphing_wing/mechanics_20261018T043837155`.
- Tet10: `results/morphing_wing/mechanics_20261018T043947293`.

Los identificadores de carpeta usan UTC; no son fechas locales de Bogotá.
Cada corrida guarda sus fuentes exactas, caso, sensores, modos e interfaces.

## Ensayo no lineal pequeño, desplazamiento −10 µm

El control total Lagrangiano SVK se ejecutó sobre el Tet4 completo con
continuación, usando fuerzas, tangentes y energía por bloques contrastados
con la implementación independiente elemento a elemento. Se ajustó la
normalización de equilibrio para un comando impuesto con fuerza exterior
nula: usa también la reacción de la frontera accionada, evitando dividir por
una carga exterior ficticia o por un valor casi cero.

La solución converge en 16 evaluaciones Newton registradas (cuatro etapas),
con mínimo det(F)=0,998512, máximo de deformación de Green 0,004752 y TE
máximo 41,9689 µm. La reacción x del servo es −1,64489 N y el máximo de von
Mises **de Cauchy** es 5,00185 MPa. La diferencia de TE respecto a diez veces
la prueba lineal Tet4 es 0,00343%, coherente con la amplitud pequeña.

Es **verificación de la ruta no lineal en el CAD nuevo**, no precisión física
aceptada del Tet4 ni resistencia comprobada. La discrepancia Tet4/Tet10 de
rigidez sigue aplicando y no se oculta con la convergencia de Newton.
La ley SVK no incluye plasticidad, daño, viscoelasticidad o adhesión entre capas.
No se extrapolan estos valores a −12,5/−25 mm.

Fuente: `results/morphing_wing/nonlinear_20261018T044418373`.

## Estado del estudio principal y de los complementos

Se actualizaron condiciones, referencias y tabla de **231 configuraciones
planificadas**, con aire de Bogotá, incidencias y fracciones del recorrido.
La tabla no incluye CL/CD/Cm inventados. Dos columnas muestran los comandos
que resultarían de elegir −12,5 o −25 mm; ninguno se declara final confirmado.
Las referencias de fuerzas y momento sobre el eje se verifican analíticamente.

El CFD 3D estacionario con transición a Reynolds de túnel no está ejecutado.
No se usan el NS laminar de bajo Re o el UVLM potencial como sustitutos de
esa validación. Deben resolverse la geometría de actuación finita, la
convergencia mecánica, el modelo de transición y la resolución/dominio fluido
antes de emitir predicciones cuantitativas del ensayo.

Las capacidades de vibración, flutter/estabilidad y resistencia se conservan.
Para el CAD nuevo hay modos y tensiones preliminares; todavía no hay umbral
de flutter ni resistencia calificada. Las campañas del STEP anterior siguen
siendo históricas y no se reinterpretan como resultados de esta semiala.

## Controles de software y revisión

La suite completa pasa **75 controles**: 69 existentes y seis nuevos sobre
referencias de flujo, momentos respecto al eje, control por desplazamiento,
selección física de apoyos/servo, sólido no lineal por bloques y una barra
analítica con alargamiento finito impuesto sin fuerza exterior.

El manifiesto explícito está en `results/morphing_wing/inputs.json`.
`run_morphing_review` genera tablas, comparación gráfica, hashes y reporte
local en `results/morphing_wing/review/`; no selecciona automáticamente la
última carpeta ni hereda el cierre de la geometría anterior.
