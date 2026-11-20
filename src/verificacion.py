import csv
from decimal import Decimal, getcontext
from pathlib import Path
import numpy as np
from gauss import (cargar_datos_texto, cargar_datos, gaussiano, sustitucion_regresiva, resolver_con_lu,)

getcontext().prec = 50

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTADOS_DIR = BASE_DIR / "resultados"

def datos_a_decimal(filas_a, filas_b):
    A_dec = [[Decimal(val) for val in fila] for fila in filas_a]
    b_dec = [Decimal(val) for val in filas_b]
    return A_dec, b_dec

#Resuelve A*x = b en Decimal para tener una referencia exacta, a diferencia de float64
def resolver_exacto_decimal(A_dec, b_dec):
    n = len(b_dec)
    M = []
    for i in range(n):
        fila = list(A_dec[i])
        fila.append(b_dec[i])
        M.append(fila)

    for k in range(n - 1):
        #Buscar pivoteo parcial en columna k
        p = k
        max_val = abs(M[k][k])
        for r in range(k + 1, n):
            if abs(M[r][k]) > max_val:
                max_val = abs(M[r][k])
                p = r

        if p != k:
            M[k], M[p] = M[p], M[k]

        for i in range(k + 1, n):
            m = M[i][k] / M[k][k]
            for j in range(k + 1, n + 1):
                M[i][j] = M[i][j] - m * M[k][j]

    x_exact = [Decimal(0)] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n]
        for j in range(i + 1, n):
            s = s - M[i][j] * x_exact[j]
        x_exact[i] = s / M[i][i]

    return x_exact

def calcular_residuo_extendido(A_dec, b_dec, x):
    x_dec = [Decimal(str(float(val))) for val in x] #Toma el valor binario exacto del float64
    n = len(b_dec)
    r_dec = []
    for i in range(n):
        ax_i = sum(A_dec[i][j] * x_dec[j] for j in range(n))
        r_dec.append(b_dec[i] - ax_i)
    norma_inf = max(abs(ri) for ri in r_dec)
    return r_dec, norma_inf

#Cálculo de error relativo real
def error_relativo_exacto(x, x_exact):
    max_diff = Decimal(0)
    max_ref = Decimal(0)
    for i in range(len(x_exact)):
        x_dec = Decimal(str(float(x[i])))
        diff = abs(x_dec - x_exact[i])
        if diff > max_diff:
            max_diff = diff
        if abs(x_exact[i]) > max_ref:
            max_ref = abs(x_exact[i])
    return float(max_diff / max_ref)

def contar_digitos_exactos(x, x_exact):
    digitos = []
    for i in range(len(x_exact)):
        x_dec = Decimal(str(float(x[i])))
        err_rel = abs(x_dec - x_exact[i]) / abs(x_exact[i])
        if err_rel == 0:
            digitos.append(15)
        else:
            d = int(np.floor(-np.log10(float(err_rel))))
            digitos.append(max(1, min(15, d)))
    return digitos

def verificar_dependencia_f7(A):
    f3 = A[2]
    f5 = A[4]
    f7 = A[6]

    # Ecuaciones normales para encontrar c1*f3 + c2*f5 que es aprox f7
    a11 = np.dot(f3, f3)
    a12 = np.dot(f3, f5)
    a22 = np.dot(f5, f5)
    b1 = np.dot(f3, f7)
    b2 = np.dot(f5, f7)

    det = a11 * a22 - a12 * a12
    c1 = (b1 * a22 - b2 * a12) / det
    c2 = (a11 * b2 - a12 * b1) / det

    f7_aprox = c1 * f3 + c2 * f5
    dist_rel = float(np.max(np.abs(f7 - f7_aprox)) / np.max(np.abs(f7)))
    return [float(c1), float(c2)], dist_rel

def refinar_solucion(A_dec, b_dec, L, U, perm, x0):
    r0_dec, norm_r0 = calcular_residuo_extendido(A_dec, b_dec, x0)
    r0_float = np.array([float(ri) for ri in r0_dec], dtype=np.float64)
    dx = resolver_con_lu(L, U, perm, r0_float)
    x1 = x0 + dx
    r1_dec, norm_r1 = calcular_residuo_extendido(A_dec, b_dec, x1)
    return x1, r0_dec, norm_r0, r1_dec, norm_r1

def analizar_sistema(A, b, x_con, x_sin, x_ref, x_exact):
    eps = float(np.finfo(np.float64).eps)
    kappa_2 = float(np.linalg.cond(A, 2))
    kappa_inf = float(np.linalg.cond(A, np.inf))
    cota_error = kappa_2 * eps
    digitos_cota = max(0.0, -float(np.log10(cota_error)))

    x_np = np.linalg.solve(A, b)
    err_rel_con = error_relativo_exacto(x_con, x_exact)
    err_rel_sin = error_relativo_exacto(x_sin, x_exact)
    err_rel_ref = error_relativo_exacto(x_ref, x_exact)
    err_rel_np = error_relativo_exacto(x_np, x_exact)
    coef_f7, dist_rel_f7 = verificar_dependencia_f7(A)

    return {"eps": eps, "kappa_2": kappa_2, "kappa_inf": kappa_inf, "cota_error": cota_error, "digitos_cota": digitos_cota, "x_np": x_np,
        "err_rel_con": err_rel_con, "err_rel_sin": err_rel_sin, "err_rel_ref": err_rel_ref, "err_rel_np": err_rel_np, "coef_f7": coef_f7, "dist_rel_f7": dist_rel_f7,}


def guardar_solucion_csv(res):
    RESULTADOS_DIR.mkdir(parents=True, exist_ok=True)
    ruta = RESULTADOS_DIR / "solucion.csv"
    x_con = res["x_con"]
    x_ref = res["x_ref"]
    dig_con = res["dig_con"]
    dig_ref = res["dig_ref"]

    with open(ruta, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["incognita", "valor_pivoteo", "digitos_pivoteo", "valor_refinado", "digitos_refinado",])
        for i in range(len(x_ref)):
            val_con = f"{x_con[i]:.{dig_con[i]}g}"
            val_ref = f"{x_ref[i]:.{dig_ref[i]}g}"
            writer.writerow([f"x{i + 1}", val_con, dig_con[i], val_ref, dig_ref[i]])


def guardar_residuo_txt(res):
    RESULTADOS_DIR.mkdir(parents=True, exist_ok=True)
    ruta = RESULTADOS_DIR / "residuo.txt"
    an = res["analisis"]

    lineas = [
        "VERIFICACION INDEPENDIENTE Y CONDICIONAMIENTO",
        "=============================================",
        "",
        "1. Condicionamiento y cota teorica (IEEE 754 float64)",
        f"Epsilon de maquina (eps): {an['eps']:.6e}",
        f"Numero de condicion kappa_2(A): {an['kappa_2']:.6e}",
        f"Numero de condicion kappa_inf(A): {an['kappa_inf']:.6e}",
        f"Cota pesimista kappa_2(A) * eps: {an['cota_error']:.6e}",
        f"Digitos garantizados por cota: {an['digitos_cota']:.2f} digitos",
        "",
        "2. Norma infinito del residuo r = b - Ax (Decimal 50 digitos)",
        f"||r||_inf sin pivoteo: {float(res['norm_r_sin']):.6e}",
        f"||r||_inf con pivoteo escalado (paso 0): {float(res['norm_r0']):.6e}",
        f"||r||_inf tras 1 refinamiento iterativo: {float(res['norm_r1']):.6e}",
        "",
        "3. Contraste de control con numpy.linalg.solve",
        f"Error relativo sin pivoteo: {an['err_rel_sin']:.6e}",
        f"Error relativo con pivoteo escalado: {an['err_rel_con']:.6e}",
        f"Error relativo tras refinamiento: {an['err_rel_ref']:.6e}",
        "",
        "4. Vector residuo r = b - Ax (Con pivoteo vs Refinado)",
    ]

    for i, (r0_i, r1_i) in enumerate(zip(res["r0"], res["r1"]), 1):
        lineas.append(
            f"r[{i:2d}] inicial = {float(r0_i):.6e} | refinado = {float(r1_i):.6e}"
        )

    ruta.write_text("\n".join(lineas) + "\n", encoding="utf-8")

def mostrar_verificacion(res):
    an = res["analisis"]
    print("VERIFICACION Y CONDICIONAMIENTO")
    print("kappa_2(A):", an["kappa_2"])
    print("Cota kappa*eps:", an["cota_error"], "-> digitos cota:", round(an["digitos_cota"], 2))
    print("||r||_inf sin pivoteo:", float(res["norm_r_sin"]))
    print("||r||_inf con pivoteo:", float(res["norm_r0"]))
    print("||r||_inf refinado:", float(res["norm_r1"]))
    print("Error rel. vs numpy (sin, con, ref):", an["err_rel_sin"], an["err_rel_con"], an["err_rel_ref"])
    print()
    for i in range(len(res["x_ref"])):
        val_con = f"{res['x_con'][i]:.{res['dig_con'][i]}g}"
        val_ref = f"{res['x_ref'][i]:.{res['dig_ref'][i]}g}"
        print(f"x{i + 1}: piv={val_con} ({res['dig_con'][i]} dig) | ref={val_ref} ({res['dig_ref'][i]} dig)")

def ejecutar_verificacion(A, b, filas_a, filas_b, U_con, info_con, x_con, x_sin):
    A_dec, b_dec = datos_a_decimal(filas_a, filas_b)
    x_exact = resolver_exacto_decimal(A_dec, b_dec)

    r_sin_dec, norm_r_sin = calcular_residuo_extendido(A_dec, b_dec, x_sin)
    x_ref, r0_dec, norm_r0, r1_dec, norm_r1 = refinar_solucion(
        A_dec, b_dec, info_con["L"], U_con, info_con["perm"], x_con
    )
    analisis = analizar_sistema(A, b, x_con, x_sin, x_ref, x_exact)
    dig_con = contar_digitos_exactos(x_con, x_exact)
    dig_ref = contar_digitos_exactos(x_ref, x_exact)

    res = {
        "x_exact": x_exact,
        "x_con": x_con,
        "x_sin": x_sin,
        "x_ref": x_ref,
        "dig_con": dig_con,
        "dig_ref": dig_ref,
        "r_sin": r_sin_dec,
        "norm_r_sin": norm_r_sin,
        "r0": r0_dec,
        "norm_r0": norm_r0,
        "r1": r1_dec,
        "norm_r1": norm_r1,
        "analisis": analisis,
    }
    guardar_solucion_csv(res)
    guardar_residuo_txt(res)
    return res

if __name__ == "__main__":
    filas_a, filas_b = cargar_datos_texto()
    A, b = cargar_datos()

    U_con, b_con, info_con = gaussiano(A, b, pivoteo_escalado=True)
    x_con, _ = sustitucion_regresiva(U_con, b_con)

    U_sin, b_sin, _ = gaussiano(A, b, pivoteo_escalado=False)
    x_sin, _ = sustitucion_regresiva(U_sin, b_sin)

    res = ejecutar_verificacion(A, b, filas_a, filas_b, U_con, info_con, x_con, x_sin)
    mostrar_verificacion(res)