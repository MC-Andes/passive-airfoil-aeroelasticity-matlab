# Datos de terceros

El código propio se distribuye bajo MIT. El extracto
`data/benchmarks/huang2026_modal.csv` y sus metadatos proceden de:

Ke Huang y Jiaying Zhang, *Experimental and numerical data from comparative
analysis of candidate structural models for active camber morphing wings*,
Mendeley Data, versión 2, [DOI 10.17632/8jxw5533py.2](https://doi.org/10.17632/8jxw5533py.2).

Licencia del dataset: [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/).
La extracción selecciona 27 observaciones de nueve geometrías, conserva las
estadísticas publicadas y añade IDs, media de tres réplicas y banderas de
auditoría. Transformación y procedencia se documentan en
`tools/import_huang2026.py` y `data/benchmarks/huang2026_provenance.json`.
Los autores originales no han avalado este solver o su uso de los datos.

Se debe conservar atribución y CC BY 4.0 al redistribuir el extracto.
La licencia MIT del código no sustituye la de estos datos.
