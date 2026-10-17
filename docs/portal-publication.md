# Publicación del estudio de aeroelasticidad en MC-Andes

La ficha se prepara en el repositorio
[del portal](https://github.com/MC-Andes/mc-andes.github.io), bajo la ruta estable
`/projects/passive-airfoil-aeroelasticity/`. Su publicación requiere revisión e
integración del PR del portal; una compilación local o un PR verde no confirma
un despliegue de GitHub Pages.

## Alcance editorial y trazabilidad

La ficha resume el estudio condicional del CAD y enlaza una revisión exacta de
este código. Distingue las 69 verificaciones del software, las comparaciones
entre discretizaciones y los ensayos de placa de terceros. No atribuye
validación del espécimen a esos controles. Las propiedades de Tough PLA,
amortiguamiento, apoyo y guía siguen siendo hipótesis autorizadas.

Las figuras son redibujos de resultados propios con la paleta del portal.
La comparación de placas deriva de Huang y Zhang, Mendeley Data v2,
[DOI 10.17632/8jxw5533py.2](https://data.mendeley.com/datasets/8jxw5533py/2),
CC BY 4.0: conserva esa licencia y la atribución de los datos. No se copian
figuras editoriales de artículos. Las otras figuras numéricas originales y el
video se publican bajo MIT. El STEP y los archivos MAT/CSV voluminosos de las
campañas permanecen fuera de Git.

En `docs/assets/airfoil/` del portal se conservan datos reducidos de publicación,
un registro de las fuentes con SHA256, unidades, tiempos del video y hashes de
los recursos. Son un extracto curado; no equivalen al archivo completo de
simulación. La licencia y advertencias de los datos importados siguen en
[THIRD_PARTY_NOTICES](../THIRD_PARTY_NOTICES.md).

## Regenerar el material visual

Primero deben existir las campañas documentadas en
[resultados del CAD](cad-original-results.md) y su manifiesto explícito
`results/cad_original/final_inputs.json`. No se eligen carpetas por su fecha ni
se sustituyen resultados faltantes con datos sintéticos. Desde la raíz:

```matlab
addpath('matlab');
export_portal_fields('/ruta/al/portal/docs/assets/airfoil');
```

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r tools/requirements-portal.txt
.venv/bin/python tools/export_portal_assets.py /ruta/al/portal/docs/assets/airfoil
```

El exportador requiere NumPy, Matplotlib, h5py, Pillow y FFmpeg distribuido por
imageio-ffmpeg; son dependencias opcionales de publicación. No intervienen en
el núcleo FEM/fluido. El video contiene 241 estados reales, no interpolados,
en η=0,5: el desplazamiento respecto al equilibrio está en micrómetros.
La ventana física 0–1,2 s se reproduce a 30 cuadros/s en 8,03 s, sin audio.
Los ejes muestran desplazamiento, no una geometría amplificada.

La geometría proviene de la tapa del sólido original de 587 caras; la estela corresponde
al promedio laminar 2D Re100 t*=12,4–14. Su vorticidad se calcula mediante
las mismas diferencias centradas usadas en el posprocesador MATLAB y se
satura explícitamente en ±8 U/c. La máscara sólida oculta campos internos.

## Comprobar y mantener la ficha

En el clon del portal: `mkdocs build --strict` y
`python tools/check_airfoil_page.py`. Verificar además escritorio/móvil,
reproducción del MP4, anclas, navegación por teclado, recursos y metadatos.
Registrar fecha, revisión del proyecto, revisión del portal y evidencias de CI.
Comprobar la URL pública después del despliegue; no afirmar publicación a
partir del artefacto de compilación.

Al actualizar cifras, regenerar figuras y hashes juntos; mantener separado el
error experimental de placa del error de discretización. Revisar permisos y
licencias antes de añadir material de terceros. Confirmar con el equipo el
responsable individual y los integrantes; la atribución corporativa MC-Andes
no permite deducir esos roles de los autores de commits.

La validación futura requiere material/proceso real, montaje/slider y datos de
presión, forma o velocidad con incertidumbre. No alterar los datos ni el
material para hacer coincidir un resultado. Cambiar hipótesis o dominio fluido
exige una nueva campaña identificada, no solo cambiar el texto de la ficha.
