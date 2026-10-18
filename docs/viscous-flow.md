# Flujo viscoso propio y alcance de la estela

El paquete MATLAB `nsflow` integra Navier–Stokes incompresible laminar en 2D/3D:
volúmenes finitos MAC, flujos convectivos conservativos centrados, difusión,
RK2 y proyección de presión. En fluido uniforme la proyección usa FFT;
con penalización resuelve el operador variable −div(beta grad p) mediante
PCG de media nula y precondicionador FFT, con beta=1/(1+dt*chi/eta).
Así la imposición de pared y la incomprensibilidad se resuelven conjuntamente.
La presión del dominio tiene media periódica nula.
Se verifican proyección/masa, Poiseuille, Taylor–Green y Beltrami ABC 3D.
El caso externo del cilindro contrasta fuerzas publicadas y sensibilidad del
dominio; no se utiliza para ajustar una constante de arrastre.

La frontera inmersa usa penalización Brinkman implícita. Es una aproximación
a no deslizamiento: ancho de máscara, tiempo de penalización, paso temporal,
resolución y dominio deben estudiarse. La caja es periódica con una franja
reservorio que repone la corriente. El balance global incluye por separado
la fuerza del reservorio y la reacción de pared. La reacción se integra como rho*chi*(u−u_pared)/eta y cierra el
balance global de cantidad de movimiento, sin calibración empírica. La variante
anterior de proyección separada se conserva para diagnósticos; las campañas
autoritativas usan `penaltyProjection='coupled'`.

`run_cad_viscous` lee la envolvente exterior auditada de la sección central.
El cuerpo ocupa el interior de esa envolvente: las cavidades se consideran
selladas para flujo exterior. En 3D se extruye con la envergadura original y
extremos tapados. Esta hipótesis conserva el exterior del álabe; no predice
presión interna ni comunicación de huecos con extremos abiertos. La incidencia
positiva del flujo equivale a girar el cuerpo un ángulo negativo en la corriente
horizontal. El sentido y la transformación se guardan en el alcance del caso.

Las magnitudes son adimensionales con U=c=rho=1 y nu=1/Re. Para interpretarlas
en aire se usa U=Re*mu/(rho*c); no se confunde Re=20/100 con U=10 m/s del caso
UVLM. El caso de heave prescrito verifica una pared móvil y lleva máscara y
velocidad de pared coherentes. No demuestra excitación espontánea de la cola.

Las salidas incluyen fuerza, presión, velocidad media/RMS, déficit y vorticidad
calculados. `surface_traction` muestrea presión y esfuerzo viscoso a una distancia
declarada fuera de la capa penalizada. Guarda Cp y tracciones, además de su
balance frente a la fuerza inmersa. Para Cp se resta la presión media del
reservorio aguas arriba y exterior al cuerpo; esa referencia se guarda.
El cambio de presión constante conserva fuerzas/momentos sobre la piel cerrada. La distancia de muestreo es una fuente de
error y se debe refinar; no se presenta ese Cp como presión de túnel validada.
La transferencia a estructura debe mantener el mismo operador cinemático y su
adjunto, con fuerzas por longitud en 2D y área en 3D.

```matlab
addpath('matlab')
run_viscous_cylinder();
run_cad_viscous(struct('Re',100,'count',[512,256],'duration',10));
run_cad_viscous(struct('Re',100,'count',[512,256], ...
    'duration',10,'amplitude',.002,'omega',2));
run_cad_viscous(struct('Re',20,'count',[192,96,96], ...
    'bounds',[-1,3;-1,1;-1.2,1.2],'duration',2));
```

El último es un transitorio 3D finito; el tiempo y los extremos deben resolverse
antes de interpretar estadísticas estacionarias. Una comparación entre dos
mallas 3D no elimina por sí sola el error de dominio o de penalización.

`cad.resource_viscous` estima un presupuesto explícito de malla uniforme para
el Reynolds del caso de referencia. Diez celdas a través de una estimación de
capa límite laminar ya superan ampliamente la RAM local. No se promete DNS a
Reynolds de túnel. URANS/LES y transición necesitan modelos y validación que
este núcleo laminar no incorpora. El UVLM sigue sirviendo para circulación y
realimentación en flujo adherido; el solver viscoso sirve para sus casos de
Reynolds/resolución declarados.

Referencias numéricas y supuestos están en las
[fichas trazables](../data/benchmarks/cad_literature.json) y el
[informe del CAD](cad-original-results.md). No hay PIV/anemometría disponibles
del espécimen; las hipótesis autorizadas no sustituyen esas mediciones.


Las campañas desarrolladas permiten reiniciar con `restartFolder` en la misma
malla. Se guarda el hash del estado inicial y la ventana estadística empieza
en `t_inicial+0.6*(t_final−t_inicial)`. Reiniciar una malla distinta se rechaza.
La comparación de las tres mallas 2D usa la ventana común 12,4–14 tiempos c/U;
los campos y perfiles proceden de medias calculadas, sin estelas dibujadas.
El [informe](cad-original-results.md) registra la dispersión de fuerzas además
de la convergencia de perfiles: una no sustituye a la otra.


`nsflow.plot_cross_wake` muestra en 3D un plano aguas abajo con déficit,
velocidad transversal y vorticidad en corriente calculados; guarda la posición
real del plano y los arrays. `run_cad_viscous_postprocess(folders,outputFolder)`
permite conservar por separado comparaciones 2D y 3D, sin mezclar Reynolds ni
sobrescribir una revisión con la otra.
