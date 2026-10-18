# Changelog

## Unreleased — geometría actualizada, PETG y servo

- Nuevo contrato de estudio con STEP/IPT identificados por hash, eje cuadrado,
  apoyo adherido y accionamiento por desplazamiento; conserva las campañas
  y capacidades anteriores como referencia/complementos.
- PETG, aire de Bogotá y campaña de 231 configuraciones planificadas, sin
  atribuir resultados a un CFD con transición todavía no implementado.
- Nuevos observadores, convenciones de CL/CD/Cm respecto al eje y pruebas
  completas del sólido actualizado.
- Sólido SVK Tet4 por bloques verificado y convergencia no lineal compatible
  con desplazamiento impuesto sin fuerza exterior; seis controles nuevos
  (75 en total).

## Unreleased — CAD original, mecánica y flujo acoplado

- Preparación de ficha científica MC-Andes, exportador de estados y figuras
  reales, guía editorial y licencias diferenciadas de los recursos.

- Mallas prismáticas del cap original, sin borrar detalles; tres niveles Tet10
  y comparación con nodos intermedios proyectados al CAD.
- Inversión Tri6 y presión sobre una huella física curva, con conservación de
  fuerza, momento y trabajo; sensores materiales sobre las pieles exteriores.
- Acción elemental estable, ordenamiento disperso, residuos de energía y
  comparación de dinámica FEM completa/reducida.
- Guías A/B/C, multiplicadores, contacto afín, fricción incremental, gravedad
  y sólido total Lagrangiano con tangente consistente verificada.
- UVLM propio, Kutta/Kelvin, estela congelada/libre y masa añadida;
  comparaciones Theodorsen/Sears, movimiento CAD prescrito y FSI fuerte.
- Eliminación exacta de circulación ligada y convolución por bloques para
  reducir coste sin truncar la historia emitida.
- Navier–Stokes MAC 2D/3D, proyección, frontera Brinkman acoplada a presión,
  pared móvil y presupuestos explícitos para evitar extrapolación a DNS de túnel.
- Hipótesis autorizadas por el usuario, fichas de literatura, bandas
  condicionales y fuentes exactas guardadas al iniciar campañas nuevas.

## Unreleased — borde de salida, 2026-10-17

- Placa MITC4 propia y cargas de presión sobre parches físicos no coincidentes.
- Interfaz modal común desde matrices de sólido/placa; presión temporal, FRF y Newmark.
- Medidas TE en coordenadas físicas, RMS, espectro y comparación sincronizada.
- Nueve geometrías PLA y 27 observaciones publicadas, con licencia y auditoría de estadísticas.
- Diez controles nuevos; ejemplos de comparación modal y vibración con presión impuesta.
- Documentada la separación entre presión prescrita, acoplamiento futuro y estela no implementada.

## Publicación inicial, 2026-10-14

- FEM sólido 3D propio con Tet4, masa consistente, apoyos y tensiones.
- Análisis modal y Newmark estructural; VLM estacionario de herraduras.
- Transferencia afín conservativa y ejecución desacoplada de una vía.
- Modelo reducido cuasiestacionario para experimentos numéricos.
- Mallas y secciones de referencia necesarias para reproducir los ejemplos.
- Referencias de Tough PLA y ABS, documentación técnica y datos experimentales requeridos.
- Verificación local y automatizada; manifiesto de integridad de entradas.
- Eliminada la dependencia implícita de `range` a toolboxes adicionales;
  el cálculo de envergadura usa máximo menos mínimo en MATLAB base.
- CI fijada en Ubuntu 24.04 para evitar cambios implícitos de plataforma.

### Limitaciones históricas de la publicación inicial

- Diferencia de 66,1% en desplazamiento máximo entre dos mallas del caso STEP.
- Material, empotramiento y tratamiento unido del mecanismo son supuestos.
- Sin contacto, grandes deformaciones FEM, CFD viscoso ni flutter validado.

## Desarrollo CAD original — 2026-10-17

- Contratos físicos y campañas trazables sobre el STEP original, sin
  sobrescribir referencias; atlas, intervalos de material y sensores de cola.
- FEM Tet10, Tri6, masa consistente, superficies curvas y transferencia
  adjunta; presión sobre un parche físico recortado y verificación de exposición.
- Equilibrado de sistemas dispersos y refinamiento del equilibrio; once
  controles adicionales, 44 en total.
- Preprocesamiento rechaza mallas degeneradas o no conformes; comprobación
  del repositorio admite archivos nuevos sin añadirlos al índice Git.
- La convergencia estructural y validación del mecanismo permanecen sujetas
  a los resultados y criterios registrados en el plan, sin anticipar FSI.
