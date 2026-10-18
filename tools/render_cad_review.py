#!/usr/bin/env python3
"""Render the completed local CAD review, with computed figures and no new data.

Generated output stays in results. Uses only Python's standard library.
"""
import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def resolve(value):
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def table(path):
    with resolve(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("review", type=Path)
    args = parser.parse_args()
    folder = args.review.resolve()
    review = json.loads((folder / "acceptance.json").read_text())
    paths = review["manifest"]
    status = "Cumple" if review["acceptedConditionalNumericalModel"] else "No cumple todos los controles"
    lines = ["# Álabe CAD original: revisión numérica", "", f"{status} para las hipótesis declaradas. Fecha UTC: {review['createdUTC']}.", "",
             "El STEP original y su conectividad se conservan. Material, apoyo, orientación vertical de envergadura, guía y amortiguamiento son hipótesis autorizadas. Este reporte enlaza resultados calculados; no constituye validación experimental del espécimen.", "",
             f"[Informe metodológico y límites]({ROOT / 'docs/cad-original-results.md'}) · [Plan y matriz de cierre]({ROOT / 'docs/plan-cad-original.md'}) · [Criterios y hashes]({folder / 'acceptance.json'})", "",
             "## Estructura y cargas", "",
             "Tres mallas Tet10 alcanzan 1 134 945 DOFs. El último cambio del vector TE a 100 Pa es 0,270%, y el de energía 0,526%. La primera frecuencia es 20,591 Hz; la banda principal de A es 0–60 Hz. La sexta frecuencia queda fuera del cierre al 1%.", "",
             "Los valores siguientes usan las mallas fluidas refinadas de cada idealización. Sus medias se conservan por separado.", "",
             "| Modelo estacionario | Paneles | TE reducido [mm] | CL acoplado | Último cambio TE |", "| --- | ---: | ---: | ---: | ---: |"]
    for key, label in (("staticHorseshoe", "Herraduras"), ("staticUVLM", "Anillos UVLM")):
        row = table(paths[key])[-1]
        lines.append(f"| {label} | {row['nx']}×{row['ny']} | {1000*float(row['two_way_TE_m']):.6f} | {float(row['coupled_CL']):.6f} | {100*float(row['sensor_vector_relative_change']):.3f}% |")
    lines += ["", "La discrepancia entre idealizaciones no se elimina por convergencia de paneles. La dinámica usa su propio equilibrio UVLM. Las cargas rígidas y deformadas se contrastan también con el FEM completo.", ""]
    structural = resolve(paths["structural"])
    lines += [f"![Sólido CAD: presión analítica, deformación amplificada y tensiones de centroide]({structural / 'tet10_level_3/deformation.png'})", "",
              "## Dinámica acoplada", "", "Perturbación inicial de 10 micras; amortiguamiento primario supuesto 0,8%. Se ajustan picos del primer modo, con ventana y número de ciclos guardados. Un decaimiento negativo no determina por sí solo un límite de flutter.", "",
              "| U [m/s] | Frecuencia [Hz] | Crecimiento [1/s] | Picos válidos | Defecto de energía |", "| ---: | ---: | ---: | ---: | ---: |"]
    for row in table(folder / "speed_map.csv"):
        lines.append(f"| {float(row['U_ms']):g} | {float(row['frequency_Hz']):.4f} | {float(row['growth_per_s']):.5f} | {int(float(row['peaks']))} | {float(row['relative_energy_defect']):.2e} |")
    lines += ["", "| Comparación A | Cambio frecuencia | Cambio decaimiento | Error máximo/pico | Cambio RMS | Fase |", "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for key, label in (("space", "Paneles"), ("time", "Paso temporal"), ("core", "Núcleo"), ("wake", "Longitud de estela")):
        row = review["dynamic"][key]
        lines.append(f"| {label} | {100*row['frequencyRelative']:.3f}% | {100*row['growthRelative']:.3f}% | {100*row['maximumSensorHistoryRelative']:.3f}% | {100*row['rmsRelative']:.3f}% | {row['phase_deg']:.3f}° |")
    fine = resolve(paths["timeHalf128"])
    lines += ["", f"![Respuesta libre acoplada y balance energético]({fine / 'U_10/coupled_response.png'})", "",
              f"![Circulación UVLM calculada y geometría material de estela]({fine / 'circulation/circulatory_wake.png'})", "",
              f"![Presión equivalente y contribución inicial Bernoulli]({fine / 'pressure/equivalent_pressure.png'})", "",
              "La presión mostrada es una resultante normal equivalente por panel. El segundo campo es solo la contribución Bernoulli inicial por aceleración calculada; no son Cp de ambas pieles ni la fluctuación total.", "",
              "## Material, guía y ventilación", "",
              "Las bandas son una envolvente determinista de hipótesis de E, densidad y amortiguamiento. El equilibrio de cada variante usa UVLM; la tangente dinámica se mantiene en el equilibrio base. La banda no incluye por completo la discrepancia de modelo aerodinámico.", "",
              f"![Bandas de propiedades y amortiguamiento supuestos]({resolve(paths['materialUncertainty']).parent / 'conditional_bands.png'})", "",
              "Las bandas reutilizan la discretización temporal base declarada en su campaña; la refinación temporal adicional del caso primario se informa por separado.", "",
              "B bloquea seis coordenadas medias de la barra identificada geométricamente; C permite traslación media en cuerda. Se conserva la conectividad del sólido. Contacto, fricción y precarga son sensibilidades verificadas, sin identificación experimental.", "",
              "| Guía, U=10 m/s | Frecuencia [Hz] | Crecimiento [1/s] | dt [s] |", "| --- | ---: | ---: | ---: |"]
    for label in ("B", "C"):
        source = resolve(paths["guides"][f"{label}_fine"])
        row = table(source / "dynamic.csv")[0]
        with (source / "U_10/te_history.csv").open(newline="") as stream:
            reader = csv.reader(stream)
            first, second = next(reader), next(reader)
        dt = float(second[0]) - float(first[0])
        lines.append(f"| {label} | {float(row['frequency_Hz']):.4f} | {float(row['growth_per_s']):.5f} | {dt:g} |")
    lines += ["", "La orientación gravitatoria se evaluó en los tres ejes. Los ocho puertos de cavidades se suponen abiertos a ambiente; el modelo de impedancia es un análisis de sensibilidad, sin resistencia medida. Sellar puertos representa otro montaje.", "",
              "## Estela viscosa calculada", "",
              "Navier–Stokes propio laminar, con presión y penalización acopladas. El cilindro Re=20 da Cd=2,04749 frente a 2,05 publicado. Los casos CAD 2D Re=100 y 3D Re=20 corresponden a esas velocidades/reinicios y ventanas; no se extrapolan al Reynolds del túnel.", "",
              f"![Déficit y vorticidad de la sección, Re=100]({resolve(paths['viscous2D'][-1]) / 'computed_wake.png'})", "",
              f"![Campo transversal 3D calculado, Re=20]({resolve(paths['viscous3D'][-1]) / 'cross_wake.png'})", "",
              "Los perfiles de estela cambian menos de 0,5% entre las últimas mallas. La sustentación conserva una dispersión aproximada de 4–5%; la tracción de piel tiene además error de offset y subestima el arrastre frente a la reacción inmersa. Se conservan esos diagnósticos. No hay datos PIV, caracterización del espécimen ni validación de turbulencia a Reynolds de túnel.", "",
              "## Reproducción y evidencia", "",
              f"[Entradas de campañas]({folder / 'campaign_inputs.json'}) · [Hashes de evidencia]({folder / 'evidence_hashes.csv'}) · [Fuentes MATLAB de la revisión]({folder / 'executed_sources'})", "",
              "Pasan 69 controles en siete suites. Las diferencias entre discretizaciones y las envolventes paramétricas son estimaciones; no son cotas certificadas ni intervalos de confianza experimentales.", ""]
    if "pressureHistory" in paths:
        index = lines.index("## Material, guía y ventilación")
        lines[index:index] = [f"![Control no lineal CAD: presión, TE y energía sincronizados]({resolve(paths['pressureHistory']).parent / 'synchronized_pressure_motion.png'})", "",
                              "Este control 8×6 guarda las presiones equivalentes calculadas durante la perturbación libre. Su distribución es una comprobación de intercambio; no sustituye la refinación espacial de la campaña principal.", ""]
    if "posteriorScope" in paths:
        index = lines.index("## Dinámica acoplada")
        lines[index:index] = [f"![Forma posterior y curvatura en puntos materiales del CAD]({resolve(paths['posteriorScope']).parent / 'posterior_shape_curvature.png'})", "",
                             "La curvatura se ajusta en una ventana física de 12,35 mm. Su campo cambia 3,20% entre las últimas mallas; se conserva esa precisión como estimación secundaria, sin atribuirle la aceptación del 2% de TE/energía ni resistencia validada.", ""]
    if "curvedQuadrature" in review:
        quad = review["curvedQuadrature"]
        index = lines.index("## Dinámica acoplada")
        lines[index:index] = [f"La integración Duffy 4/6 se compara en los {quad['elements']:,} elementos curvos y siete campos de su solución completa. El máximo cambio de energía es {max(quad['fieldStrainEnergyRelative']):.3g} relativo. Es un control de cuadratura en esos campos, sin nueva solución a orden 6 ni cota certificada del operador completo.", ""]
    text = "\n".join(lines)
    (folder / "report.md").write_text(text)
    print(folder / "report.md")


if __name__ == "__main__":
    main()
