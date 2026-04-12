/**
 * Implementação da decomposição de Cholesky para resolução de sistemas lineares
 * Ax = b, onde A é simétrica e definida positiva.
 * 
 * A decomposição de Cholesky fatora a matriz A no produto L * L^T,
 * onde L é uma matriz triangular inferior.
 * 
 * Complexidade: O(n³) para a decomposição, O(n²) para cada resolução.
 * 
 */
public class Cholesky {

    private int n;          // ordem da matriz (número de linhas/colunas)
    private double[][] L;   // fator L da decomposição (triangular inferior)

    /**
     * Exceção lançada quando a matriz não é definida positiva.
     */
    public static class NaoDefinidaPositiva extends RuntimeException {
        public NaoDefinidaPositiva(String mensagem) {
            super(mensagem);
        }
    }

    /**
     * Constrói a decomposição de Cholesky a partir da matriz A.
     * 
     * @param a matriz simétrica e definida positiva (n x n)
     * @throws IllegalArgumentException se a matriz não for quadrada
     * @throws NaoDefinidaPositiva se a matriz não for definida positiva
     */
    public Cholesky(double[][] a) {
        // Verifica se a matriz é quadrada
        if (a.length == 0) {
            throw new IllegalArgumentException("Matriz não pode ser vazia.");
        }
        for (int i = 0; i < a.length; i++) {
            if (a[i].length != a.length) {
                throw new IllegalArgumentException("Matriz deve ser quadrada.");
            }
        }

        n = a.length;
        L = new double[n][n];   // aloca o fator L (inicialmente com zeros)

        // Algoritmo de Cholesky (triangularização inferior)
        // Itera sobre cada coluna i (e linha i)
        for (int i = 0; i < n; i++) {
            // Primeiro, calcula os elementos da coluna i (da diagonal para baixo)
            for (int j = i; j < n; j++) {
                double soma = 0.0;
                // Soma sobre os elementos já calculados nas linhas anteriores (k < i)
                for (int k = 0; k < i; k++) {
                    soma += L[i][k] * L[j][k];
                }
                // Se estamos na diagonal (i == j)
                if (i == j) {
                    double valor = a[i][i] - soma;
                    if (valor <= 0.0) {
                        throw new NaoDefinidaPositiva(
                            "Matriz não é definida positiva (diagonal " + i + " <= 0)");
                    }
                    L[i][i] = Math.sqrt(valor);
                } else {
                    // Elemento fora da diagonal
                    L[j][i] = (a[j][i] - soma) / L[i][i];
                }
            }
        }
        // Nota: a matriz L é triangular inferior, e os elementos acima da diagonal
        // permanecem zero. A matriz triangular superior é a transposta de L.
    }

    /**
     * Retorna o fator L (triangular inferior) da decomposição.
     * @return cópia da matriz L (n x n)
     */
    public double[][] getL() {
        double[][] copia = new double[n][n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(L[i], 0, copia[i], 0, n);
        }
        return copia;
    }

    /**
     * Resolve o sistema A x = b utilizando os fatores L e L^T.
     * 
     * @param b vetor dos termos independentes (tamanho n)
     * @return vetor x solução (tamanho n)
     * @throws IllegalArgumentException se b tiver tamanho incorreto
     */
    public double[] resolver(double[] b) {
        if (b.length != n) {
            throw new IllegalArgumentException(
                "Vetor b tem " + b.length + " elementos, mas a matriz é " + n + "x" + n);
        }

        double[] x = new double[n];

        // --- Etapa 1: resolver L y = b  (substituição progressiva) ---
        // y é armazenado temporariamente no vetor x
        for (int i = 0; i < n; i++) {
            double soma = 0.0;
            for (int k = 0; k < i; k++) {
                soma += L[i][k] * x[k];
            }
            x[i] = (b[i] - soma) / L[i][i];
        }

        // --- Etapa 2: resolver L^T x = y  (substituição regressiva) ---
        // Agora x contém y, e vamos substituir pelo resultado final x
        for (int i = n - 1; i >= 0; i--) {
            double soma = 0.0;
            for (int k = i + 1; k < n; k++) {
                soma += L[k][i] * x[k];  // note: L^T tem elemento (i,k) = L[k][i]
            }
            x[i] = (x[i] - soma) / L[i][i];
        }
        return x;
    }

    /**
     * Calcula o produto L * y (multiplicação do fator triangular inferior por um vetor).
     * 
     * @param y vetor (tamanho n)
     * @return vetor resultante b = L y
     */
    public double[] multiplicarL(double[] y) {
        if (y.length != n) {
            throw new IllegalArgumentException("Vetor y tem tamanho incorreto.");
        }
        double[] b = new double[n];
        for (int i = 0; i < n; i++) {
            double soma = 0.0;
            for (int j = 0; j <= i; j++) {
                soma += L[i][j] * y[j];
            }
            b[i] = soma;
        }
        return b;
    }

    /**
     * Resolve o sistema L y = b, onde L é o fator triangular inferior.
     * 
     * @param b vetor dos termos independentes (tamanho n)
     * @return vetor y = L^{-1} b
     */
    public double[] resolverL(double[] b) {
        if (b.length != n) {
            throw new IllegalArgumentException("Vetor b tem tamanho incorreto.");
        }
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            double soma = 0.0;
            for (int j = 0; j < i; j++) {
                soma += L[i][j] * y[j];
            }
            y[i] = (b[i] - soma) / L[i][i];
        }
        return y;
    }

    /**
     * Calcula a inversa da matriz original A (simétrica e positiva definida).
     * 
     * @return matriz inversa (n x n)
     */
    public double[][] inversa() {
        double[][] inv = new double[n][n];
        // Resolve cada coluna da identidade para obter a inversa
        // Primeiro, resolve L * Z = I (coluna por coluna)
        // Depois resolve L^T * inv = Z (substituição regressiva)
        for (int j = 0; j < n; j++) {
            // Vetor auxiliar z (coluna da solução intermediária)
            double[] z = new double[n];
            // Resolve L * z = e_j (substituição progressiva)
            for (int i = 0; i < n; i++) {
                double soma = 0.0;
                for (int k = 0; k < i; k++) {
                    soma += L[i][k] * z[k];
                }
                double e_j = (i == j) ? 1.0 : 0.0;
                z[i] = (e_j - soma) / L[i][i];
            }
            // Agora resolve L^T * inv[:,j] = z (substituição regressiva)
            for (int i = n - 1; i >= 0; i--) {
                double soma = 0.0;
                for (int k = i + 1; k < n; k++) {
                    soma += L[k][i] * inv[j][k];  // inv[j][k] já armazena a linha j, coluna k
                }
                inv[j][i] = (z[i] - soma) / L[i][i];
                inv[i][j] = inv[j][i]; // simetria da inversa
            }
        }
        return inv;
    }

    /**
     * Calcula o logaritmo natural do determinante da matriz original A.
     * Como A = L * L^T, det(A) = (det(L))^2 = (prod_{i} L[i][i])^2.
     * Então log(det(A)) = 2 * sum(log(L[i][i])).
     * 
     * @return ln(det(A))
     */
    public double logDeterminante() {
        double soma = 0.0;
        for (int i = 0; i < n; i++) {
            soma += Math.log(L[i][i]);
        }
        return 2.0 * soma;
    }

    /**
     * Retorna o fator triangular superior R da decomposição de Cholesky,
     * tal que A = Rᵀ * R.
     * @return matriz R (n x n), triangular superior
     */
    public double[][] getR() {
        double[][] R = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                R[i][j] = L[j][i];   // L é triangular inferior
            }
        }
        return R;
    }

    public static void main(String[] args) {
        double[][] A = {
            {4, 2, 0, 0},
            {2, 5, 2, 0},
            {0, 2, 5, 2},
            {0, 0, 2 ,5}
        };
        double[] b = { 4, 6, 8, -1};

        /*
        Matriz positiva definida 4x4
        
        double[][] A = {
            {4, 2, 0, 0},
            {2, 5, 2, 0},
            {0, 2, 5, 2},
            {0, 0, 2 ,5}
        };
        double[] b = { 4, 6, 8, -1};


        */

        try {
            Cholesky cholesky = new Cholesky(A);
            
            System.out.println("Fator L (triangular inferior):");
            double[][] L = cholesky.getL();
            for (int i = 0; i < L.length; i++) {
                for (int j = 0; j < L[i].length; j++) {
                    System.out.printf("%8.4f ", L[i][j]);
                }
                System.out.println();
            }
            System.out.println("Fator R (triangular superior)");
            double[][] R = cholesky.getR();
            for (int i = 0; i < R.length; i++) {
                for (int j = 0; j < R[i].length; j++) {
                    System.out.printf("%8.4f ", R[i][j]);
                }
                System.out.println();
            }
            
            double[] x = cholesky.resolver(b);
            System.out.println("\nSolução do sistema A x = b:");
            for (int i = 0; i < x.length; i++) {
                System.out.printf("x[%d] = %.6f\n", i, x[i]);
            }
            
            System.out.printf("\nLog(det(A)) = %.6f\n", cholesky.logDeterminante());
            
            double[][] inv = cholesky.inversa();
            System.out.println("\nMatriz inversa:");
            for (int i = 0; i < inv.length; i++) {
                for (int j = 0; j < inv[i].length; j++) {
                    System.out.printf("%8.4f ", inv[i][j]);
                }
                System.out.println();
            }
            
        } catch (NaoDefinidaPositiva e) {
            System.err.println("Erro: " + e.getMessage());
        }
    }
}