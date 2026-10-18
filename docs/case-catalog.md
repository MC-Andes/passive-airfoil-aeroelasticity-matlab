# Casos publicados y comparación generalizada

Se compara por geometría, apoyos, material, carga y observable. Compartir
perfil NACA o usar PLA no hace equivalentes dos experimentos. La ejecución
automática usa placas completamente parametrizadas; los otros trabajos
orientan las extensiones al perfil y mecanismo real.

| Referencia | Relación con cola flexible | Uso actual |
| --- | --- | --- |
| [Huang y Zhang, datos v2, 2026](https://doi.org/10.17632/8jxw5533py.2) | Placas PLA con zona rígida y franja flexible posterior; ensayos modales. | Nueve geometrías y 27 observaciones incorporadas; comparación con MITC4 propio. |
| [Taguchi et al., 2020](https://doi.org/10.2322/tjsass.63.18) | Perfil pasivo con medición de deformación y aerodinámica. | Referencia para forma, cargas y presión; no se implementa una réplica de su mecanismo. |
| [Kai et al., 2023](https://doi.org/10.1007/978-981-19-2689-1_78) | NACA0012 con mecanismo compliant impreso. | Faltan entradas verificadas de piel, constitutiva, uniones y respuesta para reproducción cuantitativa. |
| [Ochi et al., 2024](https://doi.org/10.2514/1.J063000) | Simulación aeroelástica y validación de perfil pasivo impreso. | Resumen consultado; no se han extraído presiones, mallas o datos suficientes para ejecutar ese caso. |
| [Malim et al., 2023](https://doi.org/10.3390/aerospace10030256) | Palas impresas; propiedades efectivas y análisis estructural. | E y rho Tough PLA/ABS en `fem3d.material_reference`; no valida directamente vibración de cola en túnel. |

## Benchmark incorporado

Fuente: Huang y Zhang, Mendeley Data **versión 2**, DOI `10.17632/8jxw5533py.2`.
El [artículo de datos](https://doi.org/10.1016/j.dib.2026.113015) describe la
campaña. La licencia del dataset es **CC BY 4.0**, distinta de la del artículo.

El extracto conserva frecuencias adimensionales publicadas, dispersión entre
réplicas, Hz individuales y celdas originales. No incluye libros completos,
videos o modelos comerciales. URL y hash del libro quedan en
`data/benchmarks/huang2026_provenance.json`. `tools/import_huang2026.py` repite
la extracción con Python estándar y exige el hash esperado.

| AR = L/c | Fracción flexible en cuerda | Espesores nominales |
| --- | --- | --- |
| 1,8; 2,5; 5 | 0, placa uniforme | 4 mm |
| 5 | 25%; 50%; 75% | rígido 4 mm, flexible 1 mm |
| 3; 4; 6 | 50% | rígido 4 mm, flexible 1 mm |

La cuerda es 50 mm y L es longitud libre en envergadura. La franja flexible
está detrás de `x/c=1-f`; el escalón no se coloca a lo largo de la envergadura.
Se inmoviliza `y=0` de la región rígida. La raíz de la franja flexible queda
libre. Se comprobó en `AR5fef0d5PPModel.bdf`, SPC1 11, y se usa como montaje
nominal de la familia. La placa uniforme se empotra en toda su raíz.

E=2,11 GPa, nu=0,3 y rho=1180 kg/m³ proceden de MAT1 de ese archivo. Son
parámetros nominales del modelo; no propiedades medidas de todas las probetas.
Se compara `Omega=2*pi*f*L^2*sqrt(rho*tR/DR)`, con
`DR=E*tR^3/[12*(1-nu^2)]`. El campo `nominal_model_Hz` no debe compararse
directamente con Hz individuales de probetas con propiedades distintas.

## Auditoría y límites

- Test4/Test6 intercambian sus fracciones entre libro y tabla del artículo.
  Los IDs usan AR/fracción y conservan la etiqueta original.
- En filas 46–48, la media publicada usa H/L; en filas 97–99 usa D/H. La SD
  usa D/H/L. Se conservan y se añaden media de tres réplicas y bandera. No se
  corrige el libro ni se interpreta esa SD como intervalo de confianza.
- El primer modo B está identificado. Las otras comparaciones se emparejan
  por frecuencia; sin formas medidas y MAC no se confirma identidad modal.
- No se ajusta E para borrar errores. Error experimental y cambio de malla
  quedan separados. El benchmark no valida presión, estela, FSI o flutter.

En R2026a con refinaciones 8/16 en cuerda, el error del primer modo va de
+0,68% a −10,87%; el cambio entre mallas es como máximo 0,108%. La discrepancia
persiste pese al refinamiento y puede incluir propiedades, montaje, espesor e
idealización. No equivale a validar cualquier mecanismo de cola flexible.

## Agregar un caso

1. Registrar ID y fuente/versionado, ejes/unidades, región flexible, espesores,
   apoyos, slider, material y amortiguamiento. Distinguir medidas y supuestos.
2. Elegir placa, sólido u otra formulación verificada. `ael.reduce` recibe K/M
   y operadores B/O sin imponer orden de DOFs o una base modal ilustrativa.
3. Definir B sobre áreas físicas y O en los puntos medidos. Registrar área
   de cada canal de presión y estaciones del borde de salida.
4. Registrar velocidad, rho, viscosidad/Re, incidencia, paredes, excitación,
   sincronización y muestreo. Declarar escalas; no mezclar RMS y amplitud pico.
5. Comparar estática, modos y respuesta forzada. `ael.compare_series` reporta
   errores sincronizados sin extrapolar. Conservar refinaciones espacial,
   modal y temporal junto a cada caso.

La interfaz permite agregar geometrías; la evidencia disponible decide cuáles
se pueden comparar cuantitativamente. El STEP del usuario sigue pendiente de
caracterización, montaje y mecánica del slider.
