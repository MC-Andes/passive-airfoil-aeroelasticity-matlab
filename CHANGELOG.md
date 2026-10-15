# Changelog

## Unreleased — publicación inicial, 2026-10-14

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
