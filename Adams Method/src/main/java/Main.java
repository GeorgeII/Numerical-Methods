import static java.lang.Math.exp;
import static java.lang.Math.sin;
import static java.lang.Math.cos;

public class Main {

    static final int size = 3;
    static final int m = 3;
    static final double a = 0;
    static final double b = 1;
    static final double h = 0.1;
    static final int n = (int) ((b - a) / h + 1);
    static final double[][] solutionMatrixGaussianElimination = {
            {0, 2*exp(a)*sin(2*a), 2*exp(a)*cos(2*a)},
            {exp(a), -exp(a)*cos(2*a), exp(a)*sin(2*a)},
            {-exp(a), -3*exp(a)*cos(2*a), 3*exp(a)*sin(2*a)}
    };

    // x y z
    // x0 y0 z0
    // .. .. ..
    // xn yn zn
    static double[][] matrix = new double[n][];
    static double[][] funcMatrix = new double[n][];
    static double[][] preciseResult = new double[n][];
    static double[] constants;


    static double[] adams(int row, double t) {
        double[] res = new double[size];
        funcMatrix[row - 1] = func(matrix[row-1], t);
        for (int i = 0; i < size; i++) {
            res[i] = matrix[row-1][i] + h / 12 *
                    (23 * funcMatrix[row-1][i] - 16 * funcMatrix[row-2][i] + 5 * funcMatrix[row-3][i]);
        }

        return res;
    }

    static double[] func(double[] vectorOfVars, double t) {
        double[] res = new double[size];
        res[0] = vectorOfVars[0] - vectorOfVars[1] - vectorOfVars[2];
        res[1] = vectorOfVars[0] + vectorOfVars[1];
        res[2] = 3 * vectorOfVars[0] + vectorOfVars[2];

        return res;
    }

    static double[] rungeKutta(int row, double t) {
        double[][] k = new double[m][size];
        funcMatrix[row - 1] = func(matrix[row - 1], t);
        k[0] = funcMatrix[row - 1];
        double[] paramsForPass = new double[size];

        // calculate Yn + K1*h/2, where Yn is a vector (x, y, z, ...)
        for (int i = 0; i < size; i++) {
            paramsForPass[i] = matrix[row - 1][i] + h / 2 * k[0][i];
        }
        k[1] = func(paramsForPass, t + h / 2);

        // calculate Yn - K1 * h + 2*K2*h, where Yn is a vector (x, y, z, ...)
        for (int i = 0; i < size; i++) {
            paramsForPass[i] = matrix[row - 1][i] - h * k[0][i] + 2 * h * k[1][i];
        }
        k[2] = func(paramsForPass, t + h);

        double[] res = new double[size];
        for (int i = 0; i < size; i++) {
            res[i] = matrix[row - 1][i] + h / 6 * (k[0][i] + 4 * k[1][i] + k[2][i]);
        }

        return res;
    }

    static double[] preciseFunc(double t) {
        double[] res = new double[size];
        res[0] = exp(t) * (2*constants[1]*sin(2*t) + 2*constants[2]*cos(2*t));
        res[1] = exp(t) * (constants[0] - constants[1]*cos(2*t) + constants[2]*sin(2*t));
        res[2] = exp(t) * (-constants[0] - 3*constants[1]*cos(2*t) + 3*constants[2]*sin(2*t));

        return res;
    }

    public static void main(String[] args) {

        

        double t = a;

        matrix[0] = new double[size];
        // initial conditions
        matrix[0][0] = 0;
        matrix[0][1] = 1;
        matrix[0][2] = -1;
        t += h;


        constants = GaussianElimination.lsolve(solutionMatrixGaussianElimination, matrix[0]);
        for (int i = 0; i < size; i++) {
            System.out.println("Const" + (i+1) + " = " + constants[i]);
        }

        // find 1st and 2nd values
        matrix[1] = rungeKutta(1, t);
        t += h;
        matrix[2] = rungeKutta(2, t);
        t += h;

        for (int i = m; i < n; i++) {
            matrix[i] = adams(i, t);
            t += h;
        }

        t = a;
        for (int i = 0; i < n; i++) {
            preciseResult[i] = preciseFunc(t);
            t += h;
        }

        t = a;
        System.out.printf("%11s%11s%11s%11s", "t", "x", "y", "z");
        System.out.printf("%11s%11s%11s", "preciseX", "preciseY", "preciseZ");
        System.out.printf("%11s%11s%11s\n", "errorX", "errorY", "errorZ");
        for (int i = 0; i < n; i++) {
            System.out.printf("%11f", t);
            for (int j = 0; j < size; j++) {
                System.out.printf("%11.7f", matrix[i][j]);
            }
            for (int j = 0; j < size; j++) {
                System.out.printf("%11.7f", preciseResult[i][j]);
            }
            for (int j = 0; j < size; j++) {
                System.out.printf("%11.7f", Math.abs(matrix[i][j] - preciseResult[i][j]));
            }
            t += h;
            System.out.println();
        }
    }
}