# Desarrollo completado y validación posterior

La prioridad actual es el [ala morphing PETG con slider accionado](morphing-wing.md):
CAD actualizado no prismático, apoyo en eje y CFD estacionario de la semiala.
Vibraciones, flutter y resistencia se mantienen como complementos. Las campañas
anteriores pertenecen a otro STEP. El CFD 3D con transición y su convergencia
a Reynolds de túnel siguen siendo desarrollo pendiente, no una capacidad
validada por los casos laminares existentes.


El [plan del CAD original](plan-cad-original.md) se ejecuta sobre el STEP
original, con hipótesis geométricas y bibliográficas autorizadas por el usuario.
El [informe de resultados](cad-original-results.md) registra las campañas,
convergencias, errores y límites de cada modelo.

Están implementados Tet10/Tri6 curvos, transferencia física conservativa,
restricciones de guía, contacto/fricción, sólido total Lagrangiano, reducción
del FEM real, acoplamiento estático, UVLM con memoria y acoplamiento dinámico
fuerte. La extensión Navier–Stokes propia calcula estelas laminares 2D/3D con
penalización y presión acopladas. Las placas publicadas siguen siendo
comparaciones independientes del mecanismo del CAD.

La continuación experimental requiere caracterizar impresión, masa,
propiedades direccionales, amortiguamiento, montaje y funcionamiento del
slider; medir desplazamiento/presión sincronizados y estela. Esos datos
permitirían contrastar las predicciones condicionales y reducir su incertidumbre.
No constituyen entradas medidas de las campañas presentes.

Fuera del alcance verificado quedan leyes de elastómeros/viscoelasticidad sin
identificar, transición/turbulencia a Reynolds de túnel, flutter experimental y
PIV. Ampliar esos modelos requiere verificaciones y datos propios; el núcleo
laminar no se presenta como URANS, LES o DNS del túnel. El presupuesto de
recursos y las alternativas se describen en la [rama viscosa](viscous-flow.md).
