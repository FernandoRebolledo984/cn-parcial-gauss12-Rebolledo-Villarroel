import csv
from pathlib import Path
import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
DATOS_DIR = BASE_DIR / "datos"


def cargar_datos():
    ruta_a = DATOS_DIR / "matriz_a.csv"
    ruta_b = DATOS_DIR / "vector_b.csv"
    with open(ruta_a, "r", encoding="utf-8-sig") as f:
        filas_a = [[valor.strip() for valor in fila] for fila in csv.reader(f) if fila]
    with open(ruta_b, "r", encoding="utf-8-sig") as f:
        filas_b = [linea.strip() for linea in f if linea.strip()]
    A = np.array(filas_a, dtype=np.float64)
    b = np.array(filas_b, dtype=np.float64)
    return A, b


def gaussiano(A, b, pivoteo_escalado=True):
    #algoritmo de triangulizacion.
    U = A.copy()
    b_mod = b.copy()
    n = len(b_mod)
    if pivoteo_escalado:
        #pregunta si es true o false para saber si hacer el enroque de filas para pivoteo parcial escalado.
        escalas = np.max(np.abs(U), axis=1)
    for k in range(n - 1):
        if pivoteo_escalado:
            razones = np.abs(U[k:, k]) / escalas[k:]
            p = k + np.argmax(razones)
            if p != k:
                U[[k, p]] = U[[p, k]]
                b_mod[[k, p]] = b_mod[[p, k]]
                escalas[[k, p]] = escalas[[p, k]]
        for i in range(k + 1, n):
            m = U[i, k] / U[k, k]
            U[i, k] = 0.0
            U[i, k + 1:] = U[i, k + 1:] - m * U[k, k + 1:]
            b_mod[i] = b_mod[i] - m * b_mod[k]

    return U, b_mod


def sustitucion_regresiva(U, b_mod):
    n = len(b_mod)
    x = np.zeros(n, dtype=np.float64)#12 ceros
    #despeje de incognita mientras se opera con las ya conocidas desde el numero de la diagonal hacia la derecha.
    for i in range(n - 1, -1, -1):
        suma = np.dot(U[i, i + 1:], x[i + 1:])
        x[i] = (b_mod[i] - suma) / U[i, i]

    return x


if __name__ == "__main__":
    #carga datos, opera 2 veces el algoritmo sin pivoteo parcial escalado y con.(tambien hace una copia de la matriz a y vector b para no perderlos)
    A, b = cargar_datos()
    U_con, b_con = gaussiano(A, b, pivoteo_escalado=True)
    x_con = sustitucion_regresiva(U_con, b_con)
    U_sin, b_sin = gaussiano(A, b, pivoteo_escalado=False)
    x_sin = sustitucion_regresiva(U_sin, b_sin)
