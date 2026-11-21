# INFORME

## 1. Identificación

Nombres: Felipe Alonso Villarroel Ramos & Fernando Antonio Rebolledo Norambuena.  
Asignatura: Computación Numérica (Laboratorio-S2)  
Fecha: 20-11-2026  
Enlace al repositorio: https://github.com/FernandoRebolledo984/cn-parcial-gauss12-Rebolledo-Villarroel  
Uso de IA: Se utilizó la IA para mejor uso de bibliotecas (numpy, principalmente), para la lectura y salida de archivos (debido al pobre manejo de funciones estándares por nuestra parte), y también ayuda para embellecer los mensajes para los archivos de salida y usar un lenguaje más técnico. Junto con esto, se hizo uso de esta para corroborar los resultados obtenidos y, en caso de error, cambiar lo necesario.  
También se crearon las tablas presentes en este informe.

## 2. Planteamiento

Se resuelve $Ax = b$ con la matriz A (tamaño 12x12, 144 elementos no nulos y 9 cifras significativas) y b con 20 cifras significativas.
Los archivos CSV de salida se leen como texto, donde se convierten a float64 para la eliminación, y se mantienen como texto para construir "Decimal" en la verificación.
Unos de los problemas que surgieron que afectan al coste vienen siendo:
1. Hubieron operaciones de restas de múltiplos con valores absolutos muy distintos en tamaño entre sí, siendo algunos de $x \cdot 10^{5} - x \cdot 10^{-5}$ (filas 6 y 5 respectivamente), provocando que los bits bajos del operando se desplacen fuera de los 53 bits de mantisa.
2. El pivote inicial es pequeño, siendo este $a_{11} = 1.3407 \times 10^{-8}$, frente a $a_{21} = 1.19229739 \times 10^{4}$, resultando en un cociente de $8.9 \times 10^{11}$.

## 3. Estrategia de eliminación

gauss.py ejecuta tanto el pivoteo parcial escalado como sin él, siendo:
1. Con pivoteo: Elegir la fila $p \ge k$ maximizando $|a_{ik}|/s_{i}$ con $s_{i} = \max_{j}|a_{ij}|$, operando una sola vez sobre la matriz original y permutando las filas, intercambiando también las filas de $U$, de $b$, las escalas, de $L$ y del vector de permutación.
2. Para cada fila $i > k$ se calcula una vez $m_{ik} = u_{ik} / u_{kk}$, guardandolo en L, asignandod $u_{ik} = 0$ directamente y actualizar $u_{i,k + 1\dots n}$ y $b_{i}$ con el mismo $m_{ik}$.

El mínimo teórico no se excede porque se usa solo un multiplicador por fila, el barrido parte de la columna $k + 1$, existe solo una factorización sin inversa, y b se transforma en la misma pasada que A, utilizando la sustitución regresiva con un solo bucle explícito con un contador.

Es necesario igual calcular la escala una sola vez para evitar recalcularla en cada paso, evitando también cambiar el criterio de elección a medida que las filas cambian.

## 4. Conteo de operaciones

1 flop es 1 operación aritmética realizada. Debido a que las operaciones son realizadas a nivel de procesamiento, se tiene que hacer el conteo al momento de operar las matrices, realizandose en gauss.py, enviando el registro de los contadores (entre otros datos) a conteo.py. 
El sobrecosto del pivoteo es contado aparte, al igual que el cálculo inicial de las escalas $s_{i}$.

| Paso | Fila pivote (original) | Mult / Div | Sumas / Restas | Subtotal |
| :---: | :---: | :---: | :---: | :---: |
| 1 | 8 | 143 | 132 | 275 |
| 2 | 2 | 120 | 110 | 230 |
| 3 | 1 | 99 | 90 | 189 |
| 4 | 6 | 80 | 72 | 152 |
| 5 | 9 | 63 | 56 | 119 |
| 6 | 11 | 48 | 42 | 90 |
| 7 | 12 | 35 | 30 | 65 |
| 8 | 5 | 24 | 20 | 44 |
| 9 | 4 | 15 | 12 | 27 |
| 10 | 10 | 8 | 6 | 14 |
| 11 | 3 | 3 | 2 | 5 |
| **Eliminación** | | **638** | **572** | **1210** |
| **Sustitución regresiva** | | **78** | **66** | **144** |
| **Total aritmético** | | **716** | **638** | **1354** |

Las corridas con y sin pivoteo tienen la misma cantidad de flops aritméticos.

### 5. Comparación con la cota

- La cota de la pauta es 1486, el resultado es 1354, quedando en un $-8.88\%$ (dentro del $\pm 10\%$ permitido).
- La cota es calcular el elemento que se anula en cada fila (multiplicación y resta). El código lo asigna directo, siendo la diferencia $2 \times 66 = 132$, y sumado con el resultado de cota da $1354 + 132 = 1486$.

## 6. Análisis de Mantisa y Pivoteo

| Métrica | Con pivoteo escalado | Sin pivoteo |
| :--- | :---: | :---: |
| Multiplicador máximo (valor absoluto) | $2.466 \times 10^{9}$ | $2.753 \times 10^{13}$ |
| Factor de crecimiento $g$ | 1.783 | 2641.2 |
| Bits de mantisa perdidos ($\log_2 g$) | 0.83 | 11.37 |
| Error relativo final (norma $\infty$, relativo al máximo de la solución) | $2.36 \times 10^{-13}$ | $2.08 \times 10^{-10}$ |
| Máximo error relativo por componente (*) | $8.5 \times 10^{-11}$ | $7.5 \times 10^{-8}$ |
| Dígitos correctos por componente (*) | 10 a 13 | 7 a 10 |

$a_{11}$ es inaceptable como pivote porque:
1. Con $a_{11} = 1.3407 \times 10^{-8}$ los multiplicadores son enormes, siendo el mayor $m = a_{61}/a_{11} = 3.6904 \times 10^{5}/1.3407 \times 10^{-8} = 2.753 \times 10^{13}$.
2. Al restar $mF_{1}$ de la fila 6, el término $ma_{1j}$ llega a $2.55 \times 10^{9}$, mientras que la mayor entrada de la fila 6 es $9.67 \times 10^{5}$.
3. El pivoteo escalado elige en el paso 1 la fila 8 (mejor razón $|a_{i1}|/s_i$) y mantiene $g = 1.78$, perdiendo 0.83 bits en vez de $11.37$.

Entre las dos corridas la razón entre factores de crecimiento es $2641.2/1.783 = 1481$ ($\log_2 = 10.5$ bits), la razón entre errores relativos es $2.08 \times 10^{-10}/2.36 \times 10^{-13} = 882$ ($\log_2 = 9.8$ bits). 
$g$ explica la brecha de exactitud de pérdida de mantisa entre ambas corridas.

## 7. Condicionamiento

- Estimación $\kappa_2(A) = \sigma_{max}/\sigma_{min} = 2.258 \times 10^{6}/7.320 \times 10^{-8} = 3.0847 \times 10^{13}$ (no se calcula la inversa).
- Cota clásica con $\varepsilon = 2.22 \times 10^{-16}$: $\kappa_2\,\varepsilon = 6.85 \times 10^{-3}$. La cota garantiza 2.16 dígitos.
- El resultado real con pivoteo en el error relativo es $2.36 \times 10^{-13}$, con una pasada de refinamiento $6 \times 10^{-17}$

La cota resulta pesismista debido a que $\kappa_2(A)$ varía al escalado de filas; unos 9.6 de sus 13.5 órdenes de magnitud provienen de la disparidad entre escalas de fila. Al equilibrar las filas, el condificionamiento cae a $\kappa_2(D^{-1}A) = 8.1 \times 10^3$, valor que predice aproximadamente 11.7 digitos y refleja el comportamiento real del pivoteo escalado.

## 7. Resultados

| Incógnita | Valor reportado (refinado) | Dígitos afirmados | Valor con pivoteo (sin refinar) | Dígitos correctos |
| :---: | :--- | :---: | :--- | :---: |
| x1 | 0.738491562700000 | 15 | 0.7384915626999 | 13 |
| x2 | -1.20593741680000 | 15 | -1.20593741680 | 12 |
| x3 | 0.00438719520000000 | 15 | 0.004387195200 | 10 |
| x4 | -0.0917264538000000 | 15 | -0.091726453800 | 11 |
| x5 | 2.84610937250000 | 15 | 2.84610937250 | 12 |
| x6 | -0.00628154090000000 | 15 | -0.006281540899 | 10 |
| x7 | 1.47382059160000 | 15 | 1.47382059160 | 12 |
| x8 | -0.0384927164000000 | 15 | -0.038492716400 | 11 |
| x9 | 0.915264803700000 | 15 | 0.915264803700 | 12 |
| x10 | -2.60713948520000 | 15 | -2.60713948520 | 12 |
| x11 | 0.0518736294000000 | 15 | 0.051873629400 | 11 |
| x12 | -1.08493716250000 | 15 | -1.08493716250 | 12 |

Las 12 incógnitas superan los 8 dígitos correctos. Los ceros finales del valor refinado son correctos, no relleno. La solución referencia coincide con números de 10 decimales. Se afirman 15 dígitos debido a float64.

## 8. Verificación

Reutilizando L, U y la permutación sin refactorizar, una iteración reduce el error relativo en $3.9 \times 10^{3}$.

Control con biblioteca. numpy.linalg.solve (LU con pivoteo parcial) da un error de $5.23 \times 10^{-13}$, del mismo orden que la implementación propia con pivoteo ($2.36 \times 10^{-13}$, unas 2.2 veces menor). Se utilizó como control.

En resolver_exacto_decimal hace una eliminación con pivoteo parcial en Decimal de 50 dígitos, exclusivamente para verificación.

## 9. Conclusiones

Respecto al compromiso elegido, el pivoteo parcial escalado cuesta unas 143 operaciones adicionales sobre 1354 flops, pero reduce el factor de crecimiento en 1481 veces y el error en 882 veces (de $2.1 \times 10^{-10}$ a $2.4 \times 10^{-13}$), costando $O(n^{2})$ en comparación a $O(n^{3})$ de la factorización.
Con n = 1000 la eliminación costaría aproximadamente $2n^{3}$, siendo la busqueda del pivote mínima ($n^{2}/2 = 5 \times 10^{5}$ comparaciones), no se guardaría una copia de la matriz por paso y equilibraría las filas de antemano, que es barato (con $O(n^{2})$), reproduciendo lo que hace el pivoteo escalado.

## 10. Reproducibilidad

Todo el directorio resultados/ se genera ejecutando python src/main.py, con una versión de Python mínima de 3.10, numpy y módulos estándar csv, pathlib y decimal. Numpy se utiliza para álgebra básica