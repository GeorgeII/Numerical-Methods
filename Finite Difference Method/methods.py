import numpy as np
import pandas
from tabulate import tabulate


def print_tables(x, y, y_precise_solution):
    table1 = pandas.DataFrame({"x": x, "Precise y": y_precise_solution, "y": y, "Error": abs(y - y_precise_solution)})
    table1 = table1.round(6)
    print(table1)

    # tabulated output
    col_headers = ["x", "Precise y", "y", "Error"]
    merged_array = np.array([x, y_precise_solution, y, abs(y - y_precise_solution)]).T
    table2 = tabulate(merged_array , col_headers, tablefmt="fancy_grid", floatfmt = ".6f")
    print(table2)

    print("Maximum error: ", format(np.max(abs(y - y_precise_solution)), '.6f'))
