import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def main():
    x_number, t_number = 5, 4

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
    print("Boundary conditions were set: \n{}".format(y))

    a = -1 / h**2
    b = -1 / tau**2
    c = 2 / h**2 + 2 / h**2

    matrix_size = (x_number-1) * (t_number-1)
    A = np.zeros([matrix_size, matrix_size])

    # in my case function is always 0
    f = np.zeros(matrix_size)

    np.fill_diagonal(A, c)

    i, j = np.indices(A.shape)
    A[i==j] = c
    A[i==j-x_number+1] = b
    #A[i==j-1] = a
    for i in range(matrix_size-1):
        if (i+1) % (x_number-1) == 0:
            continue
        else:
            A[i][i+1] = a

    print(A)
    print(matrix_size)
    print(int((matrix_size-1) / (t_number-1)))


'''
    print('Result:\n{}'.format(y))
    plot_matrix(x, t, y)'''


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
