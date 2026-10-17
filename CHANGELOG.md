# Changelog

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

### Limitaciones conocidas

- Diferencia de 66,1% en desplazamiento máximo entre dos mallas del caso STEP.
- Material, empotramiento y tratamiento unido del mecanismo son supuestos.
- Sin contacto, grandes deformaciones FEM, CFD viscoso ni flutter validado.
