# Documentación del proyecto

El estudio principal actual calcula formas y fuerzas de una semiala morphing
accionada por servo, con CAD actualizado y PETG. Conserva vibraciones,
estabilidad/flutter y resistencia como extensiones. La
[definición actual](morphing-wing.md) sustituye hipótesis antiguas donde hay
nueva información y mantiene identificados los datos por confirmar.

## Ruta de lectura

| Necesidad | Documento |
| --- | --- |
| Estudio actual: CAD, PETG, eje, servo y campaña CFD | [Ala morphing actualizada](morphing-wing.md) |
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
| Retomar el trabajo sobre el CAD original | [Plan detallado del CAD](plan-cad-original.md) |
| Ejecutar el bloque estructural sobre el STEP original | [Campaña CAD](cad-original.md) |
| Revisar resultados y criterios del CAD | [Informe CAD](cad-original-results.md) |
| Usar mecanismos, UVLM y FSI del CAD | [Ruta acoplada](cad-coupled.md) |
| Estudiar presión, déficit y vorticidad laminares | [Flujo viscoso](viscous-flow.md) |
| Conocer trabajo pendiente | [Hoja de ruta](roadmap.md) |

## Estado científico

`run_verification` integra ocho suites: VLM/reducción, sólido, placa,
CAD/Tet10, mecánica no lineal/restricciones, UVLM/FSI, Navier–Stokes 2D/3D y seis controles nuevos del estudio morphing (75 en total).
La comparación experimental de nueve placas PLA se ejecuta por separado.
La familia prismática del CAD anterior se verifica con tres niveles y una comparación
de geometría curva; los resultados históricos no convergidos se conservan.

Los resultados del espécimen son condicionales al material, montaje, guía y
amortiguamiento supuestos. El usuario autorizó esas hipótesis mediante análisis
de la pieza y literatura; no existen ensayos de esta pieza. El UVLM no prueba
separación viscosa, y el solver laminar no se extrapola a turbulencia de túnel.
Consultar el [informe de aceptación](cad-original-results.md).
