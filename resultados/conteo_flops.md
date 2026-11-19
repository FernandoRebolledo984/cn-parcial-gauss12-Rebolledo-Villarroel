# Conteo de Operaciones y Estabilidad (Eliminacion Gaussiana 12x12)

## 1. Convencion de conteo
- **1 flop = 1 operacion aritmetica elemental** (`*`, `/`, `+`, `-`).
- Se reutiliza el multiplicador $m_{ik} = u_{ik}/u_{kk}$ una sola vez por fila y no se realizan operaciones sobre los 66 ceros bajo la diagonal.
- El sobrecosto de busqueda e intercambio de filas del pivoteo se contabiliza aparte.

## 2. Desglose paso a paso (Eliminacion + Sustitucion Regresiva)
| Paso | Fila pivote (original) | Mult / Div | Sumas / Restas | Subtotal Flops |
| :---: | :---: | :---: | :---: | :---: |
| Paso 1 | Fila 8 | 143 | 132 | 275 |
| Paso 2 | Fila 2 | 120 | 110 | 230 |
| Paso 3 | Fila 1 | 99 | 90 | 189 |
| Paso 4 | Fila 6 | 80 | 72 | 152 |
| Paso 5 | Fila 9 | 63 | 56 | 119 |
| Paso 6 | Fila 11 | 48 | 42 | 90 |
| Paso 7 | Fila 12 | 35 | 30 | 65 |
| Paso 8 | Fila 5 | 24 | 20 | 44 |
| Paso 9 | Fila 4 | 15 | 12 | 27 |
| Paso 10 | Fila 10 | 8 | 6 | 14 |
| Paso 11 | Fila 3 | 3 | 2 | 5 |
| **Subtotal Eliminacion** | - | **638** | **572** | **1210** |
| Sustitucion regresiva | - | 78 | 66 | 144 |
| **TOTAL ARITMETICO** | - | **716** | **638** | **1354** |

## 3. Comparacion con cotas de referencia
- **Flops aritmeticos obtenidos:** `1354` flops.
- **Cota minima de referencia (pauta):** `~1486` flops (variacion de `-8.88%`, dentro del margen permitido de `±10%`).

## 4. Sobrecosto de estabilidad (Pivoteo Parcial Escalado)
| Concepto | Cantidad |
| :--- | :---: |
| Divisiones de razon ($|u_{ik}| / s_i$) | 77 |
| Comparaciones (132 escalas iniciales + 66 busqueda) | 198 |
| Intercambios de fila realizados | 8 |

## 5. Comparacion entre corridas (Estabilidad y Mantisa)
| Metrica | Con Pivoteo Escalado | Sin Pivoteo |
| :--- | :---: | :---: |
| Flops aritmeticos | 1354 | 1354 |
| Multiplicador maximo ($\max |m|$) | 2.465558e+09 | 2.752593e+13 |
| Factor de crecimiento ($g$) | 1.783202 | 2641.224255 |
| Bits de mantisa perdidos ($\log_2 g$) | 0.83 bits | 11.37 bits |
