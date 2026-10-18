# Ala morphing actualizada: PETG, eje y slider accionado

El estudio principal calcula fuerzas y momentos estacionarios de una semiala
con forma impuesta por un servo. Se conservan **vibraciones, estabilidad/flutter
y resistencia mecánica como análisis complementarios**. No se eliminan sus
solvers, verificaciones o campañas anteriores. Sus resultados históricos se
mantienen ligados a su propio CAD, material y montaje.

La definición usa las nueve páginas del documento aportado por el usuario,
`Parámetros para el estudio CFD del ala morphing`, y los nuevos archivos STEP
/IPT. Las recomendaciones del documento se distinguen de datos confirmados;
los valores pendientes no se convierten en mediciones. El contrato está en
[study-data.json](../cases/morphing_wing/study-data.json); la configuración
operativa condicional, en [case.json](../cases/morphing_wing/case.json).

## Geometría y sistema de coordenadas

| Archivo / magnitud | Identificación actual |
| --- | --- |
| STEP | `Airfoil with slider alone (1).stp` |
| SHA256 STEP | `a5b289434de68e737d9c0bed552d32c642f6be22fc84d8c71d49434c501f14da` |
| SHA256 IPT | `7c729ab000e4f047765abb934bdaa4e0afb50b710b6722e180b24ad20efe1862` |
| BREP | Válido, un sólido conectado, 630 caras y 1.902 curvas |
| Volumen adaptativo | 0,0007090739570 m³, 15,1976% mayor que el CAD anterior |
| Cuerda de referencia | 0,250 m, según definición del ensayo |
| Extensión geométrica x auditada | Aproximadamente 0,246996 m; no se escala el CAD para forzar 0,250 m |
| Envergadura / área de referencia | 0,360 m / 0,090 m² |
| Alojamiento | Cuadrado de 15,2 mm; eje nominal de 15 mm y capa ideal de adhesivo de 0,1 mm por lado |

La [auditoría geométrica compacta](../data/reference/morphing_wing_geometry.json)
conserva firmas de las caras relevantes y la integral adaptativa de volumen.
El IPT nativo se identifica sin modificarlo. El cálculo utiliza el BREP del
STEP; no se afirma haber verificado la equivalencia de su árbol de operaciones
con el archivo nativo de Inventor.

El nuevo sólido **no es prismático**. Hay cambios internos a lo largo de la
envergadura, con detalles en 12/20,2/80/280/339,8/348 mm. No se utiliza el
extrusor de una sola sección de la campaña anterior. La malla nueva se genera
en 3D conservando el BREP. Las etiquetas antiguas de caras y las matrices,
modos, operadores de carga y sensores antiguos no son transferibles.

Las coordenadas canónicas son cuerda/envergadura/altura = CAD X/Z/Y.
La envergadura recorre y=0…0,360 m desde la raíz; no se centra en ±b/2 como
antes. El origen de cuerda/altura se fija a partir del centro del alojamiento
y la posición nominal c/4, z=4,3 mm del documento. Es un datum declarado; no
es un ajuste de escala o de rigidez.

## Apoyo y accionamiento utilizados

Las caras 71–74 delimitan el alojamiento cuadrado. Se restringen los tres
desplazamientos de sus nodos hasta una inserción **supuesta de 0,324 m**, la
recomendación 0,9b del documento. Las rotaciones del sólido se restringen por
el campo de desplazamientos en las cuatro caras; los nodos Tet4/Tet10 no
poseen seis grados de libertad independientes. La transición en la profundidad
de inserción depende de la discretización y se registra su última coordenada
restringida. Esa profundidad sigue sin confirmarse.

La raíz exterior y la punta quedan libres. No se conserva el empotramiento
artificial de todo x/c≤0,30 ni se pega estructuralmente la raíz al piso.
El eje y la unión adherida se idealizan rígidos; material y pared del tubo no
se dan por medidos ni se modela la capa de epóxico explícitamente.

Las caras 7–9 del STEP nuevo forman el frente coplanar de la barra inferior.
Sobre esa interfaz se impone un desplazamiento uniforme en cuerda, con
movimientos transversales nulos: es una **interfaz accionada rígida idealizada**.
No se infiere de ello un contacto distribuido validado ni se separan uniones
del sólido. Las firmas y el SHA se comprueban antes de usar esas etiquetas.
Esta elección debe contrastarse con la implementación física del servo/guía.

El estado de diseño inicial es slider cero, CAD sin deformar y sin carga
aerodinámica ni precarga prevista. Los esfuerzos residuales de impresión no
se conocen. Gravedad, contactos o compliancia adicionales deben identificarse
como variantes y no añadirse silenciosamente al estado de referencia.

## PETG y propiedades todavía condicionales

El material especificado es PETG, con relleno 100% e impresión vertical.
No se siguen usando las propiedades Tough PLA como propiedades de esta pieza.
Se toma E=1,6 GPa y ρ=1.270 kg/m³ como sustituto isotrópico de referencia,
con ν=0,35 **supuesto**, no medido. El conjunto de sensibilidad de E incluye
1,4–1,7 GPa y 3,0 GPa. El último valor procede de otra ficha que no identifica
su condición de impresión: no es una incertidumbre estadística del espécimen.

La ficha Prusament v1.1 verifica E=1,5±0,1/1,6±0,1 GPa para probetas
horizontal/vertical xz y fluencia 47±2/50±1 MPa. Sus probetas usan dos
perímetros y relleno rectilíneo; eso difiere de las paredes llenas por
perímetros propuestas para esta pieza. La ficha 3D4Makers indica E=3 GPa y
ρ=1.290 kg/m³; no se combinan esas propiedades como si procedieran de un
único material identificado.

Referencias consultadas: [Prusament PETG v1.1, ficha del fabricante alojada por distribuidor](https://www.3dprintergear.com.au/assets/files/technical-data-sheet%20(1).pdf)
y [3D4Makers PETG, ficha citada](https://cdn.shopify.com/s/files/1/0658/6459/7724/files/TDS_PETG.pdf).
La fuente metodológica del documento es Lendraitis y Lukoševičius (2025),
[DOI 10.3390/act14100498](https://doi.org/10.3390/act14100498); no se afirma
haber reproducido su ensayo.

El amortiguamiento permanece sin identificar. No se hereda ζ=0,008 del caso
PLA para presentar una predicción del PETG. Las tensiones se conservan como
complemento, pero un máximo de von Mises o la fluencia de catálogo no califican
la resistencia de paredes FDM anisótropas. Faltan ley constitutiva, variación
por orientación, daño/adhesión entre capas y convergencia regional de tensión.

Con un desplazamiento impuesto y sin carga exterior, multiplicar uniformemente
E no cambia la forma de equilibrio de un sólido homogéneo con ν fijo;
cambia la reacción y la tensión. Por tanto, no se debe inventar una banda de
formas variando solo E en ese problema. La dependencia de forma respecto a
E aparece al añadir presión/gravedad, heterogeneidad u otros cambios físicos.
Esta distinción se comprueba con una barra de solución analítica.

## Campaña aerodinámica y referencias de fuerzas

| U [m/s] | q=ρU²/2 [Pa] | Re=ρUc/μ |
| ---: | ---: | ---: |
| 15 | 99 | 184.358 |
| 25 | 275 | 307.263 |
| 35 | 539 | 430.168 |

Se usa ρ=0,88 kg/m³ y μ=1,79×10⁻⁵ Pa·s como valores de trabajo, pendientes
de las condiciones reales de ensayo. La incidencia recorre −4°…+8° en pasos
de 2°; hay once posiciones del servo. Son **231 configuraciones planificadas**,
no 231 simulaciones ya ejecutadas. El recorrido final aún no está fijado:
−12,5 y −25 mm son dos alternativas propuestas (−0,05c/−0,1c). La tabla
conserva la fracción del recorrido y las dos columnas de referencia, sin
seleccionar ninguna como límite confirmado.

CL y CD se proyectan en ejes de viento; Cm es el momento alrededor del eje
spanwise que pasa por [0,0625, 0, 0,0043] m, usando S=0,090 m² y c=0,250 m.
Con corriente [cosα,0,sinα], la dirección de sustentación es [−sinα,0,cosα].
Cm positivo alrededor de +y corresponde a cabeceo positivo. El integrador
recibe **fuerzas sobre el cuerpo en N y sus puntos físicos**, no Cp sin área.
Cambiar de referencia exige incluir el brazo de fuerza, no solo dividir un
momento por qSc.

La sección del túnel es 1×1×2 m. El dominio cercano al ensayo requiere piso
sin deslizamiento, efectos de punta y transición. AR efectiva=2,88 es una
idealización por reflexión; no reproduce la capa límite de una pared real.
La intensidad 0,6% es el objetivo/límite de diseño citado, no una medición
actual. La holgura y posición del disco siguen pendientes. Cerrar la ranura
exterior o pegar la raíz al piso para el **dominio fluido** son simplificaciones
explícitas, distintas del apoyo estructural.

El NS laminar MAC y el UVLM potencial existentes **no se reclasifican como
CFD 3D con transición validado a Re≈4,3×10⁵**. La generación de esta campaña
marca `tunnelCFDExecuted=false`. Se necesitan formulación/cierre de transición,
resolución de pared, dominio y convergencia apropiados antes de emitir polares
del túnel. No se sustituye el núcleo propio por Fluent: su mención en la
bibliografía describe otro trabajo.

## Ejecutar y comprobar los cambios

La preparación CAD requiere las dependencias opcionales ya utilizadas:
Gmsh/NumPy y cadquery-ocp para auditoría volumétrica, fuera del solver.

```bash
python tools/audit_step_volume.py "/ruta/Airfoil with slider alone (1).stp" \
  results/morphing_wing/geometry/volume_audit.json
python tools/audit_morphing_geometry.py "/ruta/Airfoil with slider alone (1).stp" \
  results/morphing_wing/geometry/signature_audit.json \
  --ipt "/ruta/Airfoil with slider alone (1).ipt"
python tools/mesh_step_volume.py "/ruta/Airfoil with slider alone (1).stp" \
  results/morphing_wing/meshes/pilot --size .006 --min-size .00015 --threads 2 \
  --volume-audit results/morphing_wing/geometry/volume_audit.json
```

```matlab
addpath('matlab');
run_morphing_study();          % define los 231 casos, no inventa resultados CFD
run_morphing_mechanics();      % prueba tangente Tet10, modos y tensiones
run_morphing_mechanics([],4);  % diagnóstico Tet4 para comparación de elementos
run_morphing_nonlinear();      % prueba pequeña SVK Tet4, no recorrido de diseño
run_verification();           % 75 controles: 69 anteriores y 6 nuevos
```

`run_morphing_mechanics` impone solo −1 µm para comprobar selección, balance,
modos y tangente. No extrapola esa solución a −12,5/−25 mm. El pequeño ensayo
SVK impone −10 µm con continuación; conserva deformación de Green, tensión
de Cauchy y determinantes. Su convergencia no valida la constitutiva PETG ni
la discretización en flexión. Las comparaciones entre elementos y mallas se
reportan en [resultados actualizados](morphing-wing-results.md).

## Complementos conservados y trabajo restante

- **Vibraciones:** modos propios del CAD nuevo con eje/slider definidos;
  FRF/transitorios requieren declarar amortiguamiento y estado de precarga.
- **Flutter/estabilidad:** se mantiene el UVLM/FSI verificado y la estimación
  de crecimiento/decaimiento. Requiere una base del CAD nuevo alrededor de
  cada equilibrio, suelo/estela coherentes y refinación independiente. Las
  tasas históricas a 2,5–10 m/s no se trasladan al ensayo nuevo de 15–35 m/s.
- **Resistencia mecánica:** se mantienen tensiones y sólido de rotación finita.
  Se deben evaluar regiones físicas, convergencia y criterios materiales
  apropiados; todavía no hay resistencia experimental ni factor de seguridad
  validado para esta impresión.

Faltan recorrido, inicio flexible en el CAD, inserción/material/pared del eje,
holgura, posición del disco, material/proceso real, ley constitutiva y
caracterización del flujo. La validación inicial recomendada es slider cero
contra mediciones de esta semiala; un CFD no se valida por coincidir con otro
CFD ni por tener un residuo pequeño.
