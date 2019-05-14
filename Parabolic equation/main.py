import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def main():
    x_number = 61
    sigma = 0.5

    # boundaries
    t_min, t_max = 0, 0.01
    x_min, x_max = 0, 1

    h = (x_max - x_min) / (x_number - 1)
    tau = sigma * h**2
    t_number = int(np.ceil((t_max - t_min) / tau)) + 1

    print('h = {}'.format(h))
    print('tau = {}'.format(tau))
    print('Matrix size: {}:{}'.format(x_number, t_number))

    x = np.linspace(x_min, x_max, x_number)
    t = np.linspace(t_min, t_max, t_number)
    y = np.zeros([x_number, t_number])
    print('x = {}'.format(x))
    print('t = {}'.format(t))

    initial_condition = (1.1 * x**2 + 2.1) * np.exp(-x)
    y[:, 0] = initial_condition

    # boundary conditions
    y[0, :] = 2.1
    y[-1, :] = 3.2 * np.exp(-1)
    np.set_printoptions(precision=2)
    print(y)

    gamma = 1
    sigma_vector = np.full(x_number - 3, sigma)
    diagonal = np.append(sigma_vector, sigma)

    for j in range(1, t_number):
        F = np.zeros(len(diagonal))
        phi = np.zeros(x_number - 2)

        for i in range(1, x_number-1):
            F[i-1] = y[i][j-1] + tau * phi[i-1]

        F[0] += sigma * y[0][j]
        F[-1] += sigma * y[-1][j]

        y[1:-1, j] = tridiagonal_matrix_algorithm(gamma * sigma_vector, -1*(1 + 2 * gamma * diagonal),
                                                           gamma * sigma_vector, -F)

    print('Result:\n{}'.format(y))
    plot_matrix(x, t, y)


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
