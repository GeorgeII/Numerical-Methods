import java.lang.Math;

public class Main {

    static final int numberOfVars = 3;
    static final double startValue = 0;
    static final double endValue = 1;
    static double t = startValue;
    static final double h = 0.1;
    static final int n = (int) ((endValue - startValue) / h);
    static final double[][] solutionMatrixGaussianElimination = {
            {0,            2*Math.exp(t)*Math.sin(2*t),   2*Math.exp(t)*Math.cos(2*t)},
            {Math.exp(t),  -Math.exp(t)*Math.cos(2*t),    Math.exp(t)*Math.sin(2*t)},
            {-Math.exp(t), -3*Math.exp(t)*Math.cos(2*t),  3*Math.exp(t)*Math.sin(2*t)}
    };




    // x  y  z
    // x0 y0 z0
    // .. .. ..
    // xn yn zn
    static double[][] matrix = new double[n][numberOfVars];
    static double[][] funcMatrix = new double[n][numberOfVars];

    static double[] constants;

    public static void main(String[] args) {

        // initial conditions
        matrix[0][0] = 0;
        matrix[0][1] = 1;
        matrix[0][2] = -1;
        t += h;

        // findConstants();
        //
        constants = GaussianElimination.lsolve(solutionMatrixGaussianElimination, matrix[0]);
        for (int i = 0; i < 3; i++) {
            System.out.println("Const" + (i+1) + " = " + constants[i]);
        }

        findFirstStep();
        findSecondStep();

        for (int i = 3; i < n; i++) {
            adams(i);
            t += h;
        }

        t = startValue;
        System.out.printf("%9s%9s%9s%9s\n", "t", "x", "y", "z");
        for (int i = 0; i < n; i++) {
            System.out.printf("%9f", t);
            for (int j = 0; j < numberOfVars; j++) {
                System.out.printf("%9.5f", matrix[i][j]);
            }
            t += h;
            System.out.println();
        }

        double[] res = new double[numberOfVars];
        /*for (int i = 0; i < numberOfVars; i++) {
            res[i] = preciseFunc(endValue, i);
            System.out.println(res[i]);
        }*/
        System.out.println("Precise x: " + preciseFunc(1, 0));
        System.out.println("Precise y: " + preciseFunc(1, 1));
        System.out.println("Precise z: " + preciseFunc(1, 2));
    }


    static void adams(int row) {
        for (int i = 0; i < numberOfVars; i++) {
            funcMatrix[row][i] = func(row - 1, i);
            matrix[row][i] = matrix[row-1][i] + h / 12 *
                    (23 * funcMatrix[row-1][i] - 16 * funcMatrix[row-2][i] + 5 * funcMatrix[row-3][i]);
        }
    }

    static double func(int row, int col) {
        // column specifies x', y', z',... function
        if (col == 0)
            return matrix[row][0] - matrix[row][1] - matrix[row][2];

        if (col == 1)
            return matrix[row][0] + matrix[row][1];

        if (col == 2)
            return 3 * matrix[row][0] + matrix[row][2];

        throw new IllegalArgumentException("Wrong parameter was passed in func");
    }

    static void findFirstStep() {
        for (int i = 0; i < numberOfVars; i++) {
            funcMatrix[0][i] = func(0, i);
            matrix[1][i] = matrix[0][i] + h * funcMatrix[0][i];
        }
        t += h;
    }

    static void findSecondStep() {
        for (int i = 0; i < numberOfVars; i++) {
            funcMatrix[1][i] = func(1, i);
            matrix[2][i] = matrix[1][i] + h * (3 / 2 * funcMatrix[1][i] - 1 / 2 * funcMatrix[0][i]);
        }
        t += h;
    }

    static void findConstants() {
        constants[0] = (3 * matrix[0][1] - matrix[0][2]) / (4 * Math.exp(startValue));

        constants[1] = (Math.cos(2*startValue) * (matrix[0][1]/Math.exp(startValue) - constants[0])/Math.sin(2*startValue) - matrix[0][0]/(2*Math.exp(startValue))) /
                (Math.cos(2*startValue) - Math.sin(2*startValue));

        constants[2] = (matrix[0][0]/(2*Math.exp(startValue)) - constants[1]*Math.sin(2*startValue)) / Math.cos(2*startValue);

        for (int i = 0; i < numberOfVars; i++) {
            System.out.println("C" + i + " = " + constants[i]);
        }
    }

    static double preciseFunc(double var, int col) {
        if (col == 0)
            //(e^x) * (2*0.196*sin(2*x) + 2*(-0.013)*cos(2*x))
            return Math.exp(var) * (2*constants[1]*Math.sin(2*var) + 2*constants[2]*Math.cos(2*var));

        if (col == 1)
            return Math.exp(var) * (constants[0] - constants[1]*Math.cos(2*var) + constants[2]*Math.sin(2*var));

        if (col == 2)
            return Math.exp(var) * (-constants[0] - 3*constants[1]*Math.cos(2*var) + 3*constants[2]*Math.sin(2*var));

        throw new IllegalArgumentException("Wrong parameter was passed in func");
    }
}