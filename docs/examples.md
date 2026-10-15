# Ejemplos, entradas y resultados

Todos los scripts resuelven sus rutas desde `mfilename('fullpath')`; el comando
`run` inicial se ejecuta desde la raíz del repositorio. Las salidas se crean en
`results/` y no forman parte del control de versiones. Se reemplazan al repetir
una corrida con el mismo script.

## 1. FEM sólido paramétrico

```matlab
addpath('matlab')
run('matlab/run_fem_demo.m')
```

| Entrada | Valor |
| --- | --- |
| Geometría | NACA00xx macizo, cuerda 0,247 m, envergadura 0,36 m, espesor relativo 0,12. |
| Discretización | 30×16×4 subdivisiones lógicas, seis tetraedros por celda. |
| Material | Tough PLA de referencia, con Poisson supuesto. |
| Apoyo | Todos los nodos con x≤0,3c empotrados, a toda la envergadura. |
| Aire | rho=1,225 kg/m³, U=10 m/s, alpha=3°. |
| VLM | 6×18 paneles sobre una superficie media plana. |
| Dinámica | Seis modos estructurales; vibración libre en el primero durante tres periodos. |

El mallado sólido trunca el perfil cerca de los extremos de cuerda para evitar
tetraedros de volumen cero. La superficie media VLM recorre la cuerda nominal.
La geometría maciza no reproduce las cavidades, flexuras ni slider de la pieza
de referencia. El empotramiento y la superficie aerodinámica son idealizaciones.

| Archivo en `results/fem_demo/` | Contenido |
| --- | --- |
| `summary.csv` | Material, tamaño de malla, masa, flujo, coeficientes, desplazamiento, tensión y primera frecuencia. |
| `modes.csv` | Frecuencias y residuales de los seis modos. |
| `energy.csv` | Tiempo y energía mecánica de la vibración libre. |
| `fem_demo.png` | Desplazamiento y tensiones, con amplificación del dibujo indicada. |
| `study.mat` | Sistema, superficie, material, solución, modos, transitorio y tabla. |

La vibración no contiene excitación ni retroalimentación del fluido. La
amplitud inicial es un dato numérico elegido para probar el integrador.

## 2. Inspección de la geometría de referencia

```matlab
run('matlab/inspect_reference.m')
```

Lee la triangulación y las envolventes incluidas. Muestra la superficie y cinco
secciones; guarda `results/reference/cad_inspection.png`. Puede abrir una
figura interactiva. La envolvente inferior presenta detalles cerca del slider
que deben revisarse antes de definir una piel aerodinámica nominal.

## 3. FEM sobre la malla de referencia

```matlab
run('matlab/run_fem_reference.m')
```

Usa 21.149 nodos y 69.082 tetraedros. Construye una superficie VLM de 4×14
paneles desde las secciones y transforma la malla al mismo origen y ejes. Las
condiciones de flujo son U=10 m/s y alpha=3°. El sólido se considera unido,
homogéneo y con el 30% anterior empotrado. No representa un contacto deslizante.

Produce `summary.csv`, `displacement_xyz_m.csv`, `fem_reference.png` y
`study.mat` bajo `results/fem_reference/`. La tabla de desplazamientos contiene
seis columnas sin encabezado: x,y,z,ux,uy,uz, todas en m.

El MAT contiene `result.meshConvergenceValidated=false` y una interpretación
preliminar. La figura indica explícitamente que la malla no está convergida.
No comparar estas tensiones con la resistencia del material como si se tratara
de una evaluación estructural validada.

## 4. Diagnóstico de convergencia del STEP

```matlab
run('matlab/run_fem_convergence.m')
```

La segunda malla tiene 85.367 nodos y 286.037 tetraedros. Se mantienen material,
superficie VLM, velocidad e incidencia. El script vuelve a construir apoyos y
transferencia de cargas para cada malla y guarda `convergence.csv`.

La corrida de referencia produjo:

| Malla | Masa [kg] | Desplazamiento máximo [mm] | Máximo von Mises [MPa] |
| --- | ---: | ---: | ---: |
| Inicial | 0,710615 | 0,00552013 | 0,373583 |
| Refinada | 0,713283 | 0,0162863 | 1,64687 |

La diferencia de desplazamiento relativa a la solución refinada es
`abs(u_fino-u_grueso)/u_fino≈0.661`. El umbral de 5% del script es un criterio
diagnóstico elegido para esta comparación, no un estándar universal. El caso
lo supera ampliamente. Las fuerzas de VLM se aplican mediante parches que
cambian con la malla, por lo que la comparación también incluye el efecto de
esa transferencia discreta sobre concentraciones locales.

## 5. Ejemplos reducidos

```matlab
run('matlab/run_demo.m')
run('matlab/run_reference.m')
```

`run_demo` ejecuta barridos, equilibrio, una linealización local y un transitorio
cuasiestacionario del modelo reducido paramétrico. Guarda CSV, PNG, FIG y MAT
en `results/demo/`.

`run_reference` usa la superficie de referencia y una estructura reducida
ilustrativa. Guarda el barrido y el estudio bajo `results/reference/`. Sus
rigideces, masas, amortiguamiento y coordenada interna del slider no se
identificaron del mecanismo. Esta ruta no constituye un reemplazo calibrado
del FEM sólido.

## Modificar un caso

Para cambiar un ejemplo se pueden editar sus entradas al principio del script.
Para un estudio propio es preferible crear una función o script independiente
que llame a la API y escriba en otro directorio. Registrar como mínimo commit,
versión de MATLAB, geometría, material, apoyos, flujo y tamaños de malla.
Consultar [API](api-reference.md) para el formato de cada estructura.
