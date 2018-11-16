import numpy as np
import pandas
import tabulate


def preciseFunc(x):
    yPreciseSolution = np.sinh(x) / np.sinh(1) - 2 * x
    print(x, yPreciseSolution)

