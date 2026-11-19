import csv
from pathlib import Path
import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
DATOS_DIR = BASE_DIR / "datos"

#==========================================================
#======================= FUNCIONES ========================
#==========================================================

#Lee datos en texto para posterior verificación.
def cargar_datos_texto():
    ruta_a = DATOS_DIR / "matriz_a.csv"
    ruta_b = DATOS_DIR / "vector_b.csv"
    with open(ruta_a, "r", encoding="utf-8-sig") as f:
        filas_a = [[valor.strip() for valor in fila] for fila in csv.reader(f) if fila]
    with open(ruta_b, "r", encoding="utf-8-sig") as f:
        filas_b = [linea.strip() for linea in f if linea.strip()]
    return filas_a, filas_b

def cargar_datos():
    filas_a, filas_b = cargar_datos_texto()
    A = np.array(filas_a, dtype=np.float64)
    b = np.array(filas_b, dtype=np.float64)
    return A, b

def gaussiano(A, b, pivoteo_escalado=True):
    #algoritmo de triangulizacion.
    U = A.copy()
    b_mod = b.copy()
    n = len(b_mod)
    L = np.zeros((n, n))
    perm = np.arange(n)
    pasos = []
    max_mult = 0.0
    max_inicial = np.max(np.abs(U))
    max_visto = max_inicial

    if pivoteo_escalado:
        #pregunta si es true o false para saber si hacer el enroque de filas para pivoteo parcial escalado.
        escalas = np.max(np.abs(U), axis=1)

    for k in range(n - 1):
        paso = {"paso":k+1, "muldiv":0, "sumres":0, "comparaciones":0, "divisiones_pivoteo":0, "intercambios":0}

        if pivoteo_escalado:
            razones = np.abs(U[k:, k]) / escalas[k:]
            paso["divisiones_pivoteo"] = razones.size #costo de estabilidad
            p = k + np.argmax(razones)
            paso["comparaciones"] = razones.size - 1
            if p != k:
                U[[k, p]] = U[[p, k]]
                b_mod[[k, p]] = b_mod[[p, k]]
                escalas[[k, p]] = escalas[[p, k]]
                L[[k, p]] = L[[p, k]]       #Multiplicador de cada fila
                perm[[k, p]] = perm[[p, k]] #Almacena la posición de cada fila original luego del intercambio del pivoteo
                paso["intercambios"] = 1

        paso["fila_pivote_original"] = int(perm[k]) + 1 #Almacena la posición original de la fila usada como pivote despues del intercambio.

        for i in range(k + 1, n):
            m = U[i, k] / U[k, k]
            L[i, k] = m #Se almacena m antes del cálculo
            U[i, k] = 0.0
            U[i, k + 1:] = U[i, k + 1:] - m * U[k, k + 1:]
            b_mod[i] = b_mod[i] - m * b_mod[k]

            ancho = U[i, k+1:].size     #Trozo de fila desde columna k+1 hasta el final
            paso["muldiv"] += 1+ancho+1 #1 division, 1 multiplicacion y el ancho de la fila.
            paso["sumres"] += ancho + 1 #Sumas/restas realizadas
            max_mult = max(max_mult, abs(m)) #Busca en multiplicador más grande con cada iteracion

        max_visto = max(max_visto, np.max(np.abs(U))) #Busca el valor más grande visto en la matriz en algun momento
        paso["matriz_aumentada"] = np.column_stack((U, b_mod))  #Matriz U con el vector b pegado a la derecha
        pasos.append(paso)

    info = {"L": L, "perm": perm, "pasos": pasos, "max_multiplicador": max_mult, "factor_crecimiento": max_visto/max_inicial}
    return U, b_mod, info

def sustitucion_regresiva(U, b_mod):
    n = len(b_mod)
    x = np.zeros(n, dtype=np.float64)   #12 ceros
    cont = {"muldiv": 0, "sumres": 0}   #Contador
    #despeje de incognita mientras se opera con las ya conocidas desde el numero de la diagonal hacia la derecha.
    for i in range(n - 1, -1, -1):
        s = b_mod[i]
        for j in range(i+1, n):
            s = s - U[i, j] * x[j]
            cont["muldiv"] += 1
            cont["sumres"] += 1
        x[i] = s / U[i, i]
        cont["muldiv"] += 1
    return x, cont

#Resuelve A x = b_nuevo reutilizando L, U y perm sin refactorizar (refinamiento en T4)
def resolver_con_lu(L, U, perm, b_nuevo):
    n = len(b_nuevo)
    c = b_nuevo[perm].astype(np.float64)    #Permutación de filas
    for i in range(1, n):   #Sustitucion progresiva con diagonal = 1
        for j in range(i):
            c[i] = c[i] - L[i, j] * c[j]
    x, x1 = sustitucion_regresiva(U, c)
    return x

def printear_cositas(nombre, info, md, sr):
    print(nombre)
    print("Max |m|: ", info["max_multiplicador"])
    print("Factor de crecimiento: ", info["factor_crecimiento"])
    print("Mult/div, sumas/restas: ", md, sr, "total", md+sr)
    return

if __name__ == "__main__":
    A, b = cargar_datos()
    for nombre, piv in (("Con piveteo escalado", True), ("SIN pivoteo", False)): #Recorre ambos casos en que se utiliza y no pivoteo.
        #U = triangular superior // b_mod = lado derecho transformado // info = diccionario con L, perm, contadores y métricas.
        U, b_mod, info = gaussiano(A, b, pivoteo_escalado=piv) 
        x, cont_sr = sustitucion_regresiva(U, b_mod)
        md = sum(p["muldiv"] for p in info["pasos"]) + cont_sr["muldiv"]
        sr = sum(p["sumres"] for p in info["pasos"]) + cont_sr["sumres"]
        printear_cositas(nombre, info, md, sr)