# Contribuir

El proyecto pertenece a [MC-Andes](https://github.com/MC-Andes). Sigue las
[normas organizacionales](https://github.com/MC-Andes/standards), con una
excepción de arranque: la importación inicial crea `main` para establecer el
código, los datos de ejecución y la CI. Las modificaciones posteriores se
proponen mediante Pull Request y revisión independiente.

## Preparar un cambio

1. Describir el problema, sus condiciones y el comportamiento esperado en un issue.
2. Crear una rama corta, por ejemplo `fix/traction-sign` o `docs/mesh-format`.
3. Mantener unidades SI y las convenciones de grados de libertad documentadas.
4. Ejecutar `python3 tools/check_repository.py` y, en MATLAB,
   `addpath('matlab'); run_verification();`.
5. Ejecutar un ejemplo relevante si el cambio afecta el flujo completo.
6. Actualizar documentación y [CHANGELOG](CHANGELOG.md).
7. Abrir un PR con Conventional Commits y explicar problema, cambio y verificación.

## Cambios numéricos

Preferir pruebas contra soluciones analíticas, invariantes físicos o datos
experimentales trazables. Un residual algebraico pequeño no demuestra
convergencia de malla. No reducir tolerancias de las pruebas para ocultar un
error ni presentar el caso STEP como convergido mientras persista la discrepancia.
Los ensayos nuevos deben identificar geometría, material, apoyos y cargas.

## Alcance del repositorio

Versionar código, documentación y entradas necesarias para ejecutar ejemplos.
Los CSV bajo `data/` son entradas numéricas, no reportes. No agregar resultados
regenerables, librerías descargadas, entornos, CAD originales ni archivos de
máquina. Si cambia una entrada, actualizar `data/manifest.json`, su procedencia
y los controles de compatibilidad. Las mallas son de desarrollo y no representan
una caracterización experimental del mecanismo.

La CI comprueba integridad y enlaces locales y ejecuta las 69 verificaciones
en MATLAB R2026a. El caso refinado completo se ejecuta localmente cuando el
cambio lo requiera; no forma parte de cada ejecución automática por su costo.

## Soporte y revisión

Usar issues del repositorio para errores reproducibles o preguntas de uso.
Incluir versión de MATLAB, sistema operativo, comando y mensaje de error.
Los autores no aprueban su propio cambio. La configuración efectiva de
protección de ramas corresponde a los mantenedores de MC-Andes.
