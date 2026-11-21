from gauss import cargar_datos_texto, cargar_datos
from conteo import (resumir_corrida, mostrar_resumen, guardar_conteo_markdown, guardar_traza_eliminacion,)
from verificacion import ejecutar_verificacion, mostrar_verificacion

def main():
    filas_a, filas_b = cargar_datos_texto()
    A, b = cargar_datos()
    res_con = resumir_corrida(A, b, pivoteo_escalado=True)
    res_sin = resumir_corrida(A, b, pivoteo_escalado=False)
    mostrar_resumen("CON PIVOTEO PARCIAL ESCALADO", res_con)
    mostrar_resumen("SIN PIVOTEO", res_sin)
    guardar_conteo_markdown(res_con, res_sin)
    guardar_traza_eliminacion(A, b, res_con, res_sin)
    res_verif = ejecutar_verificacion(A, b, filas_a, filas_b, res_con["U"], res_con["info"], res_con["x"], res_sin["x"],)
    mostrar_verificacion(res_verif)

if __name__ == "__main__":
    main()