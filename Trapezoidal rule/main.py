import numpy as np
import pandas
from tabulate import tabulate


def main():
    a = 0
    b = 2 * np.pi
    n = 10
    h = (b - a) / (n - 1)

    # trapezoidal rule
    x = np.linspace(a, b, n)
    t = np.linspace(a, b, n)
    sol_vector = (np.pi + 1) * np.cos(x)
    matrix = np.empty([n, n])
    for i in range(n):
        matrix[i] = h * np.cos(x[i] + t)

    matrix[:, 0] /= 2
    matrix[:, n-1] /= 2

    for i in range(n):
        matrix[i][i] += 1

    y = np.linalg.solve(matrix, sol_vector)
    exact_solution = np.cos(x)
    print_tables(x, y, exact_solution)
    print("n = ", n)
    print("h = ", h)

    
def print_tables(x, y, y_exact_solution):
    table1 = pandas.DataFrame({"x": x, "Precise y": y_exact_solution, "y": y, "Error": abs(y - y_exact_solution)})
    table1 = table1.round(7)
    print(table1)

    # tabulated output
    col_headers = ["x", "Exact y", "y", "Error"]
    merged_array = np.array([x, y_exact_solution, y, abs(y - y_exact_solution)]).T
    table2 = tabulate(merged_array , col_headers, tablefmt="fancy_grid", floatfmt=".7f")
    print(table2)

    print("Maximum error: ", format(np.max(abs(y - y_exact_solution)), '.7f'))


main()
