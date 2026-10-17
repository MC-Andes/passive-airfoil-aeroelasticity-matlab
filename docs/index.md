# Documentación del proyecto

El objetivo es estudiar perfiles tridimensionales fijos que se deforman
pasivamente por las cargas del aire. La implementación principal calcula las
matrices estructurales mediante un FEM sólido propio en MATLAB y aplica cargas
de un VLM de herraduras de forma desacoplada.

## Ruta de lectura

| Necesidad | Documento |
| --- | --- |
| Instalar y ejecutar desde un clon limpio | [Inicio rápido](getting-started.md) |
| Comprender componentes, unidades y grados de libertad | [Arquitectura](architecture.md) |
| Revisar la formulación estructural y sus hipótesis | [FEM sólido](solid-fem.md) |
| Consultar funciones y estructuras de datos | [Referencia de API](api-reference.md) |
| Interpretar ejemplos y archivos de salida | [Ejemplos](examples.md) |
| Preparar otra geometría y entender los CSV | [Datos y preprocesamiento](data-and-preprocessing.md) |
| Distinguir verificación, convergencia y validación | [Verificación](verification.md) |
| Usar el modelo reducido cuasiestacionario | [Modelo reducido](reduced-model.md) |
| Solicitar información al equipo experimental | [Datos experimentales](experimental-data.md) |
| Consultar bibliografía y procedencia de materiales | [Referencias](references.md) |
| Estudiar desplazamiento y vibración de la cola | [Borde de salida](trailing-edge.md) |
| Comparar geometrías y datos publicados | [Casos de referencia](case-catalog.md) |
| Conocer trabajo pendiente | [Hoja de ruta](roadmap.md) |

## Estado científico

La verificación comprende 12 controles del FEM sólido, 11 VLM/reducido y
10 de placa/dinámica/observables. La comparación experimental de nueve
placas PLA se ejecuta por separado; no valida el mecanismo del usuario.
El caso derivado del STEP conserva una diferencia de aproximadamente 66,1% en
desplazamiento máximo entre dos mallas. El volumen está auditado, pero la
discretización estructural todavía no permite una predicción cuantitativa
convergida del mecanismo real.

El material de referencia, el empotramiento anterior y el tratamiento como
sólido unido son hipótesis de desarrollo. El slider real no se identifica
automáticamente. El VLM no es CFD viscoso ni UVLM y no proporciona una velocidad
de flutter validada.
