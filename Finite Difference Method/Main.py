import numpy as np
import pandas
from tabulate import tabulate

a = 0
b = 1
h = 0.1
n = int((b - a) / h + 1)
x = np.linspace(a, b, n)
y = np.array([])
boundaryConditions = [[0, 1],
                      [0, -1]]






yPreciseSolution = np.sinh(x) / np.sinh(1) - 2 * x

# printing the table
df = pandas.DataFrame({"x": x, "yPrecise": yPreciseSolution, "z": x, "asd": yPreciseSolution})
print(df)

col_headers = ["x", "y", "z"]
merged_array = np.array([x, yPreciseSolution, x]).T
table = tabulate(merged_array , col_headers, tablefmt="fancy_grid")
print(table)
