import java.util.Arrays;

public class GramSchmidt {

    /**
     * Aplica o processo de Gram-Schmidt para converter uma base 
     * em uma base ortonormal (vetores ortogonais e de norma 1).
     * * @param A Matriz onde cada linha representa um vetor da base original.
     * @return Matriz com os vetores ortonormais.
     */
    public static double[][] ortonormalizar(double[][] A) {
        int m = A.length;    // Número de vetores
        int n = A[0].length; // Dimensão do espaço
        
        double[][] Q = new double[m][n];

        for (int i = 0; i < m; i++) {
            // 1. Começamos com o vetor original
            double[] v = A[i].clone();

            // 2. Subtraímos as projeções de v sobre os vetores já ortogonalizados
            for (int j = 0; j < i; j++) {
                double dotProduct = calcularProdutoEscalar(A[i], Q[j]);
                for (int k = 0; k < n; k++) {
                    v[k] -= dotProduct * Q[j][k];
                }
            }

            // 3. Normalizamos o vetor resultante (fazer com que tenha comprimento 1)
            double norma = calcularNorma(v);
            if (norma > 1e-15) { // Evita divisão por zero em vetores LD
                for (int k = 0; k < n; k++) {
                    Q[i][k] = v[k] / norma;
                }
            }
        }
        return Q;
    }

    private static double calcularProdutoEscalar(double[] u, double[] v) {
        double soma = 0;
        for (int i = 0; i < u.length; i++) {
            soma += u[i] * v[i];
        }
        return soma;
    }

    private static double calcularNorma(double[] v) {
        return Math.sqrt(calcularProdutoEscalar(v, v));
    }

    public static void main(String[] args) {
        // Exemplo: Dois vetores no R2
        double[][] baseOriginal = {
            {3, 1},
            {2, 2}
        };

        double[][] baseOrtonormal = ortonormalizar(baseOriginal);

        System.out.println("Base Ortonormal Resultante:");
        for (double[] vetor : baseOrtonormal) {
            System.out.println(Arrays.toString(vetor));
        }
    }
}