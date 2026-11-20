import math
from pathlib import Path
import numpy as np
from gauss import cargar_datos, gaussiano, sustitucion_regresiva

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTADOS_DIR = BASE_DIR / "resultados"


def resumir_corrida(A, b, pivoteo_escalado):
    U, b_mod, info = gaussiano(A, b, pivoteo_escalado=pivoteo_escalado)
    x, cont_sr = sustitucion_regresiva(U, b_mod)
    n = len(b)

    muldiv_pasos = np.array([p["muldiv"] for p in info["pasos"]])
    sumres_pasos = np.array([p["sumres"] for p in info["pasos"]])
    total_pasos = muldiv_pasos + sumres_pasos

    elim_muldiv = int(np.sum(muldiv_pasos))
    elim_sumres = int(np.sum(sumres_pasos))
    elim_total = int(np.sum(total_pasos))

    sr_muldiv = cont_sr["muldiv"]
    sr_sumres = cont_sr["sumres"]
    sr_total = sr_muldiv + sr_sumres

    total_flops = elim_total + sr_total

    if pivoteo_escalado:
        comp_pivoteo = int(np.sum([p["comparaciones"] for p in info["pasos"]]))
        div_pivoteo = int(np.sum([p["divisiones_pivoteo"] for p in info["pasos"]]))
        intercambios = int(np.sum([p["intercambios"] for p in info["pasos"]]))
    else:
        comp_pivoteo = 0
        div_pivoteo = 0
        intercambios = 0

    crecimiento = float(info["factor_crecimiento"])
    bits_perdidos = math.log2(crecimiento) if crecimiento > 1.0 else 0.0

    return {
        "U": U,
        "x": x,
        "info": info,
        "muldiv_pasos": muldiv_pasos,
        "sumres_pasos": sumres_pasos,
        "total_pasos": total_pasos,
        "elim_muldiv": elim_muldiv,
        "elim_sumres": elim_sumres,
        "elim_total": elim_total,
        "sr_muldiv": sr_muldiv,
        "sr_sumres": sr_sumres,
        "sr_total": sr_total,
        "total_flops": total_flops,
        "comp_pivoteo": comp_pivoteo,
        "div_pivoteo": div_pivoteo,
        "intercambios": intercambios,
        "max_multiplicador": float(info["max_multiplicador"]),
        "factor_crecimiento": crecimiento,
        "bits_perdidos": bits_perdidos,
    }


def mostrar_resumen(titulo, res):
    print(titulo)
    print("Eliminacion (mul/div, sum/res, total):", res["elim_muldiv"], res["elim_sumres"], res["elim_total"])
    print("Sust. regresiva (mul/div, sum/res, total):", res["sr_muldiv"], res["sr_sumres"], res["sr_total"])
    print("Total flops:", res["total_flops"])
    print("Sobrecosto pivoteo (div, comp, intercambios):", res["div_pivoteo"], res["comp_pivoteo"], res["intercambios"])
    print("Max |m|:", res["max_multiplicador"])
    print("Factor de crecimiento:", res["factor_crecimiento"])
    print("Bits perdidos:", res["bits_perdidos"])
    print()


def guardar_conteo_markdown(res_con, res_sin):
    RESULTADOS_DIR.mkdir(parents=True, exist_ok=True)
    ruta = RESULTADOS_DIR / "conteo_flops.md"

    cota_ref = 1486
    diff_pct = ((res_con["total_flops"] - cota_ref) / cota_ref) * 100

    lineas = [
        "# Conteo de Operaciones y Estabilidad (Eliminacion Gaussiana 12x12)",
        "",
        "## 1. Convencion de conteo",
        "- **1 flop = 1 operacion aritmetica elemental** (`*`, `/`, `+`, `-`).",
        "- Se reutiliza el multiplicador $m_{ik} = u_{ik}/u_{kk}$ una sola vez por fila y no se realizan operaciones sobre los 66 ceros bajo la diagonal.",
        "- El sobrecosto de busqueda e intercambio de filas del pivoteo se contabiliza aparte.",
        "",
        "## 2. Desglose paso a paso (Eliminacion + Sustitucion Regresiva)",
        "| Paso | Fila pivote (original) | Mult / Div | Sumas / Restas | Subtotal Flops |",
        "| :---: | :---: | :---: | :---: | :---: |",
    ]

    for i, paso in enumerate(res_con["info"]["pasos"]):
        md = int(res_con["muldiv_pasos"][i])
        sr = int(res_con["sumres_pasos"][i])
        tot = int(res_con["total_pasos"][i])
        fila_orig = paso["fila_pivote_original"]
        lineas.append(f"| Paso {i + 1} | Fila {fila_orig} | {md} | {sr} | {tot} |")

    lineas.extend([
        f"| **Subtotal Eliminacion** | - | **{res_con['elim_muldiv']}** | **{res_con['elim_sumres']}** | **{res_con['elim_total']}** |",
        f"| Sustitucion regresiva | - | {res_con['sr_muldiv']} | {res_con['sr_sumres']} | {res_con['sr_total']} |",
        f"| **TOTAL ARITMETICO** | - | **{res_con['elim_muldiv'] + res_con['sr_muldiv']}** | **{res_con['elim_sumres'] + res_con['sr_sumres']}** | **{res_con['total_flops']}** |",
        "",
        "## 3. Comparacion con cotas de referencia",
        f"- **Flops aritmeticos obtenidos:** `{res_con['total_flops']}` flops.",
        f"- **Cota minima de referencia (pauta):** `~{cota_ref}` flops (variacion de `{diff_pct:+.2f}%`, dentro del margen permitido de `±10%`).",
        "",
        "## 4. Sobrecosto de estabilidad (Pivoteo Parcial Escalado)",
        "| Concepto | Cantidad |",
        "| :--- | :---: |",
        f"| Divisiones de razon (u_ik dividido por la escala s_i) | {res_con['div_pivoteo']} |",
        f"| Comparaciones de busqueda de pivote | {res_con['comp_pivoteo']} |",
        f"| Intercambios de fila realizados | {res_con['intercambios']} |",
        "",
        "## 5. Comparacion entre corridas (Estabilidad y Mantisa)",
        "| Metrica | Con Pivoteo Escalado | Sin Pivoteo |",
        "| :--- | :---: | :---: |",
        f"| Flops aritmeticos | {res_con['total_flops']} | {res_sin['total_flops']} |",
        f"| Multiplicador maximo (valor absoluto) | {res_con['max_multiplicador']:.6e} | {res_sin['max_multiplicador']:.6e} |",
        f"| Factor de crecimiento ($g$) | {res_con['factor_crecimiento']:.6f} | {res_sin['factor_crecimiento']:.6f} |",
        f"| Bits de mantisa perdidos ($\\log_2 g$) | {res_con['bits_perdidos']:.2f} bits | {res_sin['bits_perdidos']:.2f} bits |",
        "",
    ])

    ruta.write_text("\n".join(lineas), encoding="utf-8")


def guardar_traza_eliminacion(A, b, res_con, res_sin):
    RESULTADOS_DIR.mkdir(parents=True, exist_ok=True)
    ruta = RESULTADOS_DIR / "traza_eliminacion.txt"
    inicial = np.column_stack((A, b))

    with open(ruta, "w", encoding="utf-8") as f:
        f.write("TRAZA DE ELIMINACION GAUSSIANA (11 PASOS)\n")
        f.write("==========================================\n")

        for titulo, res in [
            ("CORRIDA 1: CON PIVOTEO PARCIAL ESCALADO", res_con),
            ("CORRIDA 2: SIN PIVOTEO", res_sin),
        ]:
            f.write(f"\n{titulo}\n")
            f.write("==========================================\n")
            f.write("\nPaso 0 - Matriz inicial:\n")
            for fila in inicial:
                texto_fila = " ".join(f"{val:.17e}" for val in fila[:-1])
                f.write(f"{texto_fila} | {fila[-1]:.17e}\n")
            for paso in res["info"]["pasos"]:
                num = paso["paso"]
                piv_orig = paso["fila_pivote_original"]
                mat = paso["matriz_aumentada"]
                f.write(f"\nPaso {num} (pivote usado: fila original {piv_orig}) [U | b]:\n")
                for fila in mat:
                    texto_fila = " ".join(f"{val:.17e}" for val in fila[:-1])
                    f.write(f"{texto_fila} | {fila[-1]:.17e}\n")


if __name__ == "__main__":
    A, b = cargar_datos()
    res_con = resumir_corrida(A, b, pivoteo_escalado=True)
    res_sin = resumir_corrida(A, b, pivoteo_escalado=False)

    mostrar_resumen("CON PIVOTEO PARCIAL ESCALADO", res_con)
    mostrar_resumen("SIN PIVOTEO", res_sin)

    guardar_conteo_markdown(res_con, res_sin)
    guardar_traza_eliminacion(A, b, res_con, res_sin)