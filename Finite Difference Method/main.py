import numpy as np
import methods as ms

a = 0
b = 1
h = 0.1
n = int((b - a) / h + 1)
x = np.linspace(a, b, n)
y = np.zeros(n)
# boundary conditions
y[0] = 0
y[n-1] = -1

# tridiagonal matrix algorithm
A = np.ones(n)
C = np.ones(n)
B = np.full(n, -2 - h**2)
D = 2*x * h**2
MU = np.zeros(n)
V = np.copy(y)
for i in range(1, n-1):
    MU[i] = -C[i] / (A[i]*MU[i-1] + B[i])
    V[i] = (D[i] - A[i] * V[i - 1]) / (A[i] * MU[i - 1] + B[i])
for i in range(n-2, 0, -1):
    y[i] = MU[i]*y[i+1]+V[i]

y_precise_solution = np.sinh(x) / np.sinh(1) - 2 * x
ms.print_tables(x, y, y_precise_solution)
