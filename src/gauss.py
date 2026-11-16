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


if __name__ == "__main__":
    a, b = cargar_datos()
    print("Matriz A:", a)
    print("Vector b:", b)