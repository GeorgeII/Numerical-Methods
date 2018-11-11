import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import static java.lang.Math.cos;
import static java.lang.Math.exp;
import static java.lang.Math.sin;

public class Adams {

    private final int m = 3; // only 3-step Adams method

    private int size;
    private double a = 0;
    private double b = 1;
    private double h = 0.1;
    private int n = (int) ((b - a) / h + 1);

    private RealMatrix matrix;
    private RealMatrix funcMatrix;
    private RealMatrix preciseResult;
    private double[] constants;


    private RealMatrix solutionMatrixGaussianElimination;

    public Adams(int size) {
        this.size = size;
    }

    public void setInitialValues(double[] values) {
        matrix = new Array2DRowRealMatrix(n, size);
        matrix.setRow(0, values);
    }

    public void setSolutionMatrix(double[][] solutionMatrix) {
        solutionMatrixGaussianElimination = new Array2DRowRealMatrix(solutionMatrix);
    }

    // equations must be set manually
    private double[] func(double[] vectorOfVars, double t) {
        double[] res = new double[size];
        res[0] = vectorOfVars[0] - vectorOfVars[1] - vectorOfVars[2];
        res[1] = vectorOfVars[0] + vectorOfVars[1];
        res[2] = 3 * vectorOfVars[0] + vectorOfVars[2];

        return res;
    }

    // quations must be set manually
    private double[] preciseFunc(double t) {
        double[] res = new double[size];
        res[0] = exp(t) * (2*constants[1]*sin(2*t) + 2*constants[2]*cos(2*t));
        res[1] = exp(t) * (constants[0] - constants[1]*cos(2*t) + constants[2]*sin(2*t));
        res[2] = exp(t) * (-constants[0] - 3*constants[1]*cos(2*t) + 3*constants[2]*sin(2*t));

        return res;
    }

    private void rungeKutta() {
        double t = a + h;
        RealMatrix k = new Array2DRowRealMatrix(m, size);

        for (int i = 1; i < m; i++) {
            funcMatrix.setRow(i - 1, func(matrix.getRow(i - 1), t));
            k.setRowVector(0, funcMatrix.getRowVector(i - 1));
            RealVector paramsForPass = new ArrayRealVector(size);

            // calculate Yn + K1*h/2, where Yn is a vector (x, y, z, ...)
            paramsForPass.setSubVector(0, matrix.getRowVector(i - 1).add(k.getRowVector(0).mapMultiply(h / 2)));
            k.setRow(1, func(paramsForPass.toArray(), t + h / 2));

            // calculate Yn - K1 * h + 2*K2*h, where Yn is a vector (x, y, z, ...)
            paramsForPass.setSubVector(0, matrix.getRowVector(i - 1).subtract(k.getRowVector(0).mapMultiply(h))
                                                .add(k.getRowVector(1).mapMultiply(2 * h)));
            k.setRow(2, func(paramsForPass.toArray(), t + h));

            // write down the result vector
            matrix.setRowVector(i, matrix.getRowVector(i - 1).add((k.getRowVector(0)
                                    .add(k.getRowVector(1).mapMultiply(4)).add(k.getRowVector(2))).mapMultiply(h / 6)));

            t += h;
        }

    }

    public void solve() {
        funcMatrix = new Array2DRowRealMatrix(n, size);
        preciseResult = new Array2DRowRealMatrix(n, size);

        constants = GaussianElimination.lsolve(solutionMatrixGaussianElimination.getData(), matrix.getRow(0));
        rungeKutta();

        double t = a + h * m;
        for (int i = m; i < n; i++) {
            RealVector res = new ArrayRealVector(size);
            funcMatrix.setRow(i - 1, func(matrix.getRow(i - 1), t));
            res.setSubVector(0, matrix.getRowVector(i - 1).add(funcMatrix.getRowVector(i -1).mapMultiply(23)
                    .subtract(funcMatrix.getRowVector(i - 2).mapMultiply(16)).add(funcMatrix.getRowVector(i - 3).mapMultiply(5)).mapMultiply(h/12)));

            matrix.setRowVector(i, res);
            t += h;
        }

        t = a;
        for (int i = 0; i < n; i++) {
            preciseResult.setRow(i, preciseFunc(t));
            t += h;
        }
    }


    public void print() {
        double t = a;
        System.out.printf("%11s%11s%11s%11s", "t", "x", "y", "z");
        System.out.printf("%11s%11s%11s", "precise x", "precise y", "precise z");
        System.out.printf("%11s%11s%11s\n", "error x", "error y", "error z");
        for (int i = 0; i < n; i++) {
            System.out.printf("%11f", t);
            for (int j = 0; j < size; j++) {
                System.out.printf("%11.7f", matrix.getEntry(i, j));
            }
            for (int j = 0; j < size; j++) {
                System.out.printf("%11.7f", preciseResult.getEntry(i, j));
            }
            for (int j = 0; j < size; j++) {
                System.out.printf("%11.7f", Math.abs(matrix.getEntry(i, j) - preciseResult.getEntry(i, j)));
            }
            t += h;
            System.out.println();
        }
    }

    public void setStartAndEnd(double start, double end) {
        a = start;
        b = end;
        n = (int) ((b - a) / h + 1);
    }

    public void setStep(double step) {
        h = step;
        n = (int) ((b - a) / h + 1);
    }
}
