import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def main():
    x_number, t_number = 20,20

    # boundaries
    t_min, t_max = 0, 1
    x_min, x_max = 0, 1

    h = (x_max - x_min) / (x_number - 1)
    tau = (t_max - t_min) / (t_number - 1)

    print('h = {}'.format(h))
    print('tau = {}'.format(tau))
    print('Matrix size: {}:{}'.format(x_number, t_number))

    x = np.linspace(x_min, x_max, x_number)
    t = np.linspace(t_min, t_max, t_number)
    y = np.zeros([x_number, t_number])
    print('x = {}'.format(x))
    print('t = {}'.format(t))

    # boundary conditions
    y[0, :] = 30 * np.sin(np.pi * t)
    y[-1, :] = 20 * t
    y[:, 0] = 20 * x
    y[:, -1] = 30 * x * (1 - x)

    np.set_printoptions(precision=2)
    print("Boundary conditions are set: \n{}".format(y))

    a = -1 / h**2
    b = -1 / tau**2
    c = 2 / h**2 + 2 / tau**2

    matrix_size = (x_number-2) * (t_number-2)
    A = np.zeros([matrix_size, matrix_size])

    np.fill_diagonal(A, c)
    i, j = np.indices(A.shape)
    A[i == j-(x_number-2)] = b
    A[i == j-1] = a
    for i in range(matrix_size-1):
        if (i+1) % (x_number-2) == 0:
            A[i][i+1] = 0
    # symmetrize the matrix
    A = A + A.T - np.diag(A.diagonal())
    print(A)

    # in my case function is always 0
    f = np.zeros([x_number-2, t_number-2])
    # vector F in linear system: Ay = F
    F = np.zeros([x_number-2, t_number-2])

    for i in range(x_number-2):
        for j in range(t_number-2):
            if i == 0:
                F[i][j] -= a * y[0][j+1]
            if i == x_number - 3:
                F[i][j] -= a * y[-1][j+1]
            if j == 0:
                F[i][j] -= b * y[i+1][0]
            if j == t_number - 3:
                F[i][j] -= b * y[i+1][-1]
        F[i][j] += f[i][j]
    F = F.flatten()

    result = np.linalg.solve(A, F)
    result = np.reshape(result, (x_number-2, t_number-2))

    for i in range(1, x_number-1):
        for j in range(1, t_number-1):
            y[i][j] = result[i-1][j-1]

    print(result)
    print(y)
    plot_matrix(x, t, y)


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
