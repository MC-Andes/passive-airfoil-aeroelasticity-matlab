# Trabajo pendiente

La publicación inicial establece un solver verificable de desarrollo. Las
extensiones siguientes no se presentan como capacidades ya implementadas.

## Prioridad estructural

1. Resolver la convergencia del caso de referencia: investigar tetraedros
   cuadráticos y láminas/vigas para flexuras y pieles, con benchmarks de flexión.
2. Definir regiones físicas, espesores y orientación material de la pieza real.
3. Incorporar guías, restricciones entre cuerpos, contacto y fricción del slider
   cuando se conozca su arquitectura y modo de bloqueo.
4. Implementar grandes desplazamientos/rotaciones con fuerzas internas y
   tangente coherentes; añadir una ley material apropiada si se usan elastómeros.

## Cargas y acoplamiento

1. Controlar la huella de cargas distribuidas durante el refinamiento.
2. Validar el VLM con datos de fuerza/presión en el régimen de Reynolds de ensayo.
3. Acoplar en dos vías la estructura FEM con actualización de geometría y
   criterios de convergencia de cargas y desplazamientos.
4. Incorporar memoria aerodinámica mediante UVLM o un solver viscoso cuando
   el objetivo dinámico lo requiera. No inferir flutter de la linealización
   cuasiestacionaria actual.

## Reproducibilidad y validación

- Medir propiedades efectivas de la impresión y amortiguamiento.
- Registrar condiciones de montaje y sincronización experimental.
- Añadir comparación con benchmarks experimentales publicados, conservando
  licencias y procedencia de los datos.
- Ampliar plataformas/versiones de MATLAB verificadas.
- Empaquetar reducción modal desde el FEM con el reordenamiento e interpolación
  correctos; hoy la ruta reducida requiere una entrada estructural independiente.

Cada extensión debe aportar una comparación verificable y documentar su
dominio de validez antes de emplearse para predicciones del álabe.
