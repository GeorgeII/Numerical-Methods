import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def main():
    sigma = 0.5
    x_number = 31
    k = 3

    # boundaries
    t_min, t_max = 0, 0.2
    x_min, x_max = 0, 1

    h = (x_max - x_min) / (x_number - 1)
    tau = sigma * h
    t_number = int(np.ceil((t_max - t_min) / tau)) + 1

    print('h = {}'.format(h))
    print('tau = {}'.format(tau))
    print('Matrix size: {}:{}'.format(x_number, t_number))

    x = np.linspace(x_min, x_max, x_number)
    t = np.linspace(t_min, t_max, t_number)
    y = np.zeros([x_number, t_number])
    print('x = {}'.format(x))
    print('t = {}'.format(t))

    # in my case function is always 0
    f = np.zeros(len(x))

    initial_values = k * np.sin((k * np.pi * x) / x_max)
    initial_values_derivative = k * x * np.sin(k * np.pi * x / x_max)
    y[:, 0] = initial_values
    y[:, 1] = y[:, 0] + tau * initial_values_derivative + \
              0.5 * tau**2 * (second_order_derivative(y[:, 0], h) + f)

    # boundary conditions
    y[0, :] = 0
    y[-1, :] = 0
    np.set_printoptions(precision=2)
    print(y)

    sigma_vector = np.full(x_number - 3, sigma)
    diagonal = np.append(sigma_vector, sigma)

    for j in range(2, t_number):
        F = np.zeros(len(diagonal))
        y_difference = second_order_derivative(y[:, j-1], h)

        for i in range(1, x_number - 1):
            F[i-1] = 2 * y[i][j-1] - y[i][j-2]

        y[1:-1, j] = tridiagonal_matrix_algorithm(sigma_vector ** 2,
                                                  -1*(1 + 2 * diagonal **2),
                                                  sigma_vector ** 2, -F)

    print('Result:\n{}'.format(y))
    plot_matrix(x, t, y)


def second_order_derivative(y, h):
    result = np.zeros(len(y))
    for i in range(1, len(y) - 1):
        result[i] = (y[i+1] - 2 * y[i] + y[i+1]) / h**2

    return result


def tridiagonal_matrix_algorithm(a, b, c, d):
    nf = len(d)
    ac, bc, cc, dc = map(np.array, (a, b, c, d))
    for it in range(1, nf):
        mc = ac[it - 1] / bc[it - 1]
        bc[it] = bc[it] - mc * cc[it - 1]
        dc[it] = dc[it] - mc * dc[it - 1]

    xc = bc
    xc[-1] = dc[-1] / bc[-1]

    for il in range(nf - 2, -1, -1):
        xc[il] = (dc[il] - cc[il] * xc[il + 1]) / bc[il]

    return xc


def plot_matrix(x, t, y):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    T, X = np.meshgrid(t, x)

    surf = ax.plot_surface(X + 1e5, T + 1e5, y, cmap='autumn', cstride=2, rstride=2)
    ax.set_xlabel("x-Label")
    ax.set_ylabel("t-Label")
    ax.set_zlabel("z-Label")
    ax.set_zlim(y.min(), y.max())

    plt.show()


if __name__ == '__main__':
    main()
