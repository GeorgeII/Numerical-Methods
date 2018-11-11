import static java.lang.Math.cos;
import static java.lang.Math.exp;
import static java.lang.Math.sin;

public class NewMain {
    public static void main(String[] args) {

        int size = 3;
        double a = 0;
        double b = 1;
        double h = 0.1;
        double[] initialValues = {0, 1, -1};

        final double[][] solutionMatrixForGaussianElimination = {
                {0, 2*exp(a)*sin(2*a), 2*exp(a)*cos(2*a)},
                {exp(a), -exp(a)*cos(2*a), exp(a)*sin(2*a)},
                {-exp(a), -3*exp(a)*cos(2*a), 3*exp(a)*sin(2*a)}
        };

        Adams adms = new Adams(size);
        adms.setStartAndEnd(a, b);
        adms.setStep(h);
        adms.setInitialValues(initialValues);
        adms.setSolutionMatrix(solutionMatrixForGaussianElimination);
        adms.solve();
        adms.print();

    }
}
