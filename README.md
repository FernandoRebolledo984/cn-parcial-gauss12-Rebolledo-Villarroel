# Eliminación gaussiana 12x12 (El pivote)

Se resuelvió un sistema de 12 ecuaciones con 12 incógnitas usando eliminación gaussiana, con y sin pivoteo. Además se contaron las operaciones que hace el algoritmo y nos cercioramos que la solución sea correcta.

## Autoría

- Felipe Alonso Villarroel Ramos
- Fernando Antonio Rebolledo Norambuena

Asignatura: Computación Numérica (Laboratorio-S2), Universidad Católica del Maule.

## Requisitos

- Python 3.10 o superior
- Instalar las librerías con:

```
pip install -r requirements.txt
```

## Ejecución

Desde la carpeta principal del repositorio:

```
python src/main.py
```

Esto muestra un resumen por pantalla y crea estos archivos en `resultados/`:

- `conteo_flops.md`: conteo de operaciones paso a paso
- `traza_eliminacion.txt`: la matriz después de cada paso de la eliminación
- `solucion.csv`: las 12 incógnitas con sus dígitos correctos
- `residuo.txt`: verificación de la solución y número de condición

## Estructura de carpetas

```
cn-parcial-gauss12-Rebolledo-Villarroel/
├─ readme.md          
├─ informe.md         # informe con el análisis
├─ requirements.txt   # librerías necesarias
├─ datos/             # la matriz A y el vector b que los códigos leen
├─ src/               # los códigos (gauss, conteo, verificación y main)
└─ resultados/        # salida de datos de los códigos
```

## Resultado resumido

Solución final (con pivoteo y una pasada de refinamiento), a 10 dígitos significativos:

| Incógnita | Valor |
| :---: | :--- |
| x1 | 0.7384915627 |
| x2 | -1.205937417 |
| x3 | 0.004387195200 |
| x4 | -0.09172645380 |
| x5 | 2.846109372 |
| x6 | -0.006281540900 |
| x7 | 1.473820592 |
| x8 | -0.03849271640 |
| x9 | 0.9152648037 |
| x10 | -2.607139485 |
| x11 | 0.05187362940 |
| x12 | -1.084937162 |

- Total de operaciones (flops): 1354
- Norma infinito del residuo: 2.6e-10

## Informe

El análisis completo está en [informe.md](informe.md).

## Licencia

Trabajo de uso exclusivamente académico para el curso de Computación Numérica/Ingeniería Civil Informática/Universidad Católica del Maule.