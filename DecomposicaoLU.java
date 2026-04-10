/**
 * Algoritmos para decomposição LU de matrizes quadradas.
 * 
 * Decomposição LU: A = L * U, onde L é triangular inferior (com diagonal 1)
 * e U é triangular superior.
 * 
 * Com pivotamento parcial: PA = L * U (P é matriz de permutação de linhas)
 * Com pivotamento completo: PAQ = L * U (P e Q permutam linhas e colunas)
 * 
 * @author
 */
public class DecomposicaoLU {

    // ------------------------------------------------------------
    // 1. Decomposição LU simples (sem pivotamento)
    // ------------------------------------------------------------
    /**
     * Calcula a decomposição A = L * U (sem pivotamento).
     * A deve ser não singular.
     * 
     * @param A Matriz quadrada (n x n). Não será modificada.
     * @return Um array de duas matrizes [L, U] onde L é triangular inferior
     *         (com diagonal 1) e U é triangular superior.
     * @throws IllegalArgumentException se A não for quadrada ou for singular.
     */
    public static double[][][] decomposicaoLUSimples(double[][] A) {
        int n = A.length;
        // Verificar se é quadrada
        for (int i = 0; i < n; i++) {
            if (A[i].length != n) {
                throw new IllegalArgumentException("Matriz não é quadrada.");
            }
        }
        
        // Copiar A para U (será modificada)
        double[][] U = new double[n][n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, U[i], 0, n);
        }
        
        double[][] L = new double[n][n];
        // Inicializar L com diagonal 1 e zeros fora
        for (int i = 0; i < n; i++) {
            L[i][i] = 1.0;
        }
        
        // Algoritmo de Doolittle
        for (int k = 0; k < n - 1; k++) {
            // Verificar singularidade (pivô zero)
            if (U[k][k] == 0) {
                throw new IllegalArgumentException("Matriz singular: pivô zero na posição (" + k + "," + k + ")");
            }
            for (int i = k + 1; i < n; i++) {
                // Fator L[i][k] = U[i][k] / U[k][k]
                double fator = U[i][k] / U[k][k];
                L[i][k] = fator;
                // Atualizar linha i de U: subtrair fator * linha k
                for (int j = k; j < n; j++) {
                    U[i][j] -= fator * U[k][j];
                }
            }
        }
        // Verificar último pivô
        if (U[n-1][n-1] == 0) {
            throw new IllegalArgumentException("Matriz singular.");
        }
        
        return new double[][][]{L, U};
    }
    
    // ------------------------------------------------------------
    // 2. Decomposição PA = LU (pivotamento parcial por linhas)
    // ------------------------------------------------------------
    /**
     * Calcula a decomposição PA = L * U, onde P é matriz de permutação de linhas.
     * 
     * @param A Matriz quadrada (n x n). Não será modificada.
     * @return Um array contendo:
     *         [0] = matriz L (triangular inferior, diagonal 1)
     *         [1] = matriz U (triangular superior)
     *         [2] = vetor de permutação de linhas P (tamanho n)
     *               A permutação é tal que P * A = L * U.
     * @throws IllegalArgumentException se A não for quadrada ou for singular.
     */
    public static Object[] decomposicaoPALU(double[][] A) {
        int n = A.length;
        // Verificar quadrada
        for (int i = 0; i < n; i++) {
            if (A[i].length != n) {
                throw new IllegalArgumentException("Matriz não é quadrada.");
            }
        }
        
        // Copiar A para U (será modificada)
        double[][] U = new double[n][n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, U[i], 0, n);
        }
        
        double[][] L = new double[n][n];
        for (int i = 0; i < n; i++) {
            L[i][i] = 1.0;
        }
        
        // Vetor de permutação de linhas: inicialmente identidade
        int[] perm = new int[n];
        for (int i = 0; i < n; i++) perm[i] = i;
        
        for (int k = 0; k < n - 1; k++) {
            // --- Pivotamento parcial: encontrar linha com maior |U[i][k]|
            int maxIndex = k;
            double maxVal = Math.abs(U[k][k]);
            for (int i = k + 1; i < n; i++) {
                if (Math.abs(U[i][k]) > maxVal) {
                    maxVal = Math.abs(U[i][k]);
                    maxIndex = i;
                }
            }
            if (maxVal == 0) {
                throw new IllegalArgumentException("Matriz singular (coluna " + k + " totalmente zero)");
            }
            
            // Trocar linhas em U e em L (e no vetor perm)
            if (maxIndex != k) {
                // Trocar linhas de U
                double[] tempU = U[k];
                U[k] = U[maxIndex];
                U[maxIndex] = tempU;
                // Trocar linhas de L (apenas elementos já preenchidos)
                double[] tempL = L[k];
                L[k] = L[maxIndex];
                L[maxIndex] = tempL;
                // Trocar permutação
                int tempP = perm[k];
                perm[k] = perm[maxIndex];
                perm[maxIndex] = tempP;
            }
            
            // Se o pivô é zero, singular
            if (U[k][k] == 0) {
                throw new IllegalArgumentException("Matriz singular (pivô zero após pivotamento)");
            }
            
            // Eliminação
            for (int i = k + 1; i < n; i++) {
                double fator = U[i][k] / U[k][k];
                L[i][k] = fator;
                for (int j = k; j < n; j++) {
                    U[i][j] -= fator * U[k][j];
                }
            }
        }
        // Verificar último pivô
        if (U[n-1][n-1] == 0) {
            throw new IllegalArgumentException("Matriz singular.");
        }
        
        return new Object[]{L, U, perm};
    }
    
    // ------------------------------------------------------------
    // 3. Decomposição PAQ = LU (pivotamento completo: linhas e colunas)
    // ------------------------------------------------------------
    /**
     * Calcula a decomposição PAQ = L * U, onde P e Q são matrizes de permutação
     * de linhas e colunas, respectivamente.
     * 
     * @param A Matriz quadrada (n x n). Não será modificada.
     * @return Um array contendo:
     *         [0] = matriz L (triangular inferior, diagonal 1)
     *         [1] = matriz U (triangular superior)
     *         [2] = vetor de permutação de linhas P (tamanho n)
     *         [3] = vetor de permutação de colunas Q (tamanho n)
     *         A permutação é tal que P * A * Q = L * U.
     * @throws IllegalArgumentException se A não for quadrada ou for singular.
     */
    public static Object[] decomposicaoPAQLU(double[][] A) {
        int n = A.length;
        // Verificar quadrada
        for (int i = 0; i < n; i++) {
            if (A[i].length != n) {
                throw new IllegalArgumentException("Matriz não é quadrada.");
            }
        }
        
        // Copiar A para U (será modificada)
        double[][] U = new double[n][n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, U[i], 0, n);
        }
        
        double[][] L = new double[n][n];
        for (int i = 0; i < n; i++) {
            L[i][i] = 1.0;
        }
        
        // Vetores de permutação: inicialmente identidade
        int[] permLinhas = new int[n];
        int[] permColunas = new int[n];
        for (int i = 0; i < n; i++) {
            permLinhas[i] = i;
            permColunas[i] = i;
        }
        
        for (int k = 0; k < n - 1; k++) {
            // --- Pivotamento completo: buscar maior elemento em módulo na submatriz U[k..n-1][k..n-1]
            int maxRow = k;
            int maxCol = k;
            double maxVal = Math.abs(U[k][k]);
            for (int i = k; i < n; i++) {
                for (int j = k; j < n; j++) {
                    double val = Math.abs(U[i][j]);
                    if (val > maxVal) {
                        maxVal = val;
                        maxRow = i;
                        maxCol = j;
                    }
                }
            }
            if (maxVal == 0) {
                throw new IllegalArgumentException("Matriz singular (submatriz a partir de (" + k + "," + k + ") é nula)");
            }
            
            // Trocar linhas se necessário
            if (maxRow != k) {
                double[] tempU = U[k];
                U[k] = U[maxRow];
                U[maxRow] = tempU;
                // Trocar linhas de L (apenas os elementos já preenchidos)
                double[] tempL = L[k];
                L[k] = L[maxRow];
                L[maxRow] = tempL;
                // Atualizar permutação de linhas
                int tempP = permLinhas[k];
                permLinhas[k] = permLinhas[maxRow];
                permLinhas[maxRow] = tempP;
            }
            
            // Trocar colunas se necessário
            if (maxCol != k) {
                // Trocar colunas em U (para todas as linhas)
                for (int i = 0; i < n; i++) {
                    double temp = U[i][k];
                    U[i][k] = U[i][maxCol];
                    U[i][maxCol] = temp;
                }
                // Trocar colunas em L? Não, L é triangular inferior e só tem elementos nas colunas < i.
                // Mas precisamos também trocar colunas em L? Na verdade, as permutações de colunas afetam
                // a decomposição, mas não precisamos modificar L porque L é construída com base na matriz
                // após permutação de colunas. O algoritmo padrão para pivotamento completo aplica as trocas
                // de colunas em U, e as trocas de linhas em L e U. O vetor permColunas é atualizado.
                // L permanece consistente.
                // Atualizar permutação de colunas
                int tempQ = permColunas[k];
                permColunas[k] = permColunas[maxCol];
                permColunas[maxCol] = tempQ;
            }
            
            // Se o pivô é zero (não deveria acontecer pois maxVal>0)
            if (U[k][k] == 0) {
                throw new IllegalArgumentException("Matriz singular (pivô zero)");
            }
            
            // Eliminação
            for (int i = k + 1; i < n; i++) {
                double fator = U[i][k] / U[k][k];
                L[i][k] = fator;
                for (int j = k; j < n; j++) {
                    U[i][j] -= fator * U[k][j];
                }
            }
        }
        // Verificar último pivô
        if (U[n-1][n-1] == 0) {
            throw new IllegalArgumentException("Matriz singular.");
        }
        
        return new Object[]{L, U, permLinhas, permColunas};
    }
    
    // ------------------------------------------------------------
    // Métodos auxiliares para exibir matrizes e vetores
    // ------------------------------------------------------------
    public static void imprimirMatriz(double[][] M, String nome) {
        System.out.println(nome + ":");
        for (double[] linha : M) {
            for (double val : linha) {
                System.out.printf("%8.4f ", val);
            }
            System.out.println();
        }
    }
    
    public static void imprimirVetor(int[] v, String nome) {
        System.out.print(nome + ": [");
        for (int i = 0; i < v.length; i++) {
            System.out.print(v[i] + (i < v.length-1 ? ", " : ""));
        }
        System.out.println("]");
    }
    
    // ------------------------------------------------------------
    // Exemplo de uso e teste
    // ------------------------------------------------------------
    public static void main(String[] args) {
        // Exemplo: matriz 3x3
        double[][] A = {
            {2, 1, -1},
            {1, 2, 1},
            {1, 1, 1}
        };
        
        System.out.println("Matriz original A:");
        imprimirMatriz(A, "A");
        
        try {
            // 1. LU simples
            System.out.println("\n--- Decomposição LU simples ---");
            double[][][] lu = decomposicaoLUSimples(A);
            imprimirMatriz(lu[0], "L");
            imprimirMatriz(lu[1], "U");
            // Verificar: L*U deve ser igual a A
            double[][] LUprod = multiplicarMatrizes(lu[0], lu[1]);
            imprimirMatriz(LUprod, "L * U");
            
            // 2. PA = LU
            System.out.println("\n--- Decomposição PA = LU (pivotamento parcial) ---");
            Object[] paLu = decomposicaoPALU(A);
            double[][] Lp = (double[][]) paLu[0];
            double[][] Up = (double[][]) paLu[1];
            int[] permP = (int[]) paLu[2];
            imprimirMatriz(Lp, "L");
            imprimirMatriz(Up, "U");
            imprimirVetor(permP, "P (linhas permutadas)");
            // Verificar: P*A deve ser igual a L*U
            double[][] PA = aplicarPermutacaoLinhas(A, permP);
            double[][] LUprod2 = multiplicarMatrizes(Lp, Up);
            System.out.println("P * A:");
            imprimirMatriz(PA, "P*A");
            System.out.println("L * U:");
            imprimirMatriz(LUprod2, "L*U");
            
            // 3. PAQ = LU
            System.out.println("\n--- Decomposição PAQ = LU (pivotamento completo) ---");
            Object[] paqLu = decomposicaoPAQLU(A);
            double[][] Lc = (double[][]) paqLu[0];
            double[][] Uc = (double[][]) paqLu[1];
            int[] permPc = (int[]) paqLu[2];
            int[] permQc = (int[]) paqLu[3];
            imprimirMatriz(Lc, "L");
            imprimirMatriz(Uc, "U");
            imprimirVetor(permPc, "P (linhas)");
            imprimirVetor(permQc, "Q (colunas)");
            // Verificar: P*A*Q deve ser igual a L*U
            double[][] PAQ = aplicarPermutacaoLinhasColunas(A, permPc, permQc);
            double[][] LUprod3 = multiplicarMatrizes(Lc, Uc);
            System.out.println("P * A * Q:");
            imprimirMatriz(PAQ, "PAQ");
            System.out.println("L * U:");
            imprimirMatriz(LUprod3, "L*U");
            
        } catch (IllegalArgumentException e) {
            System.err.println("Erro: " + e.getMessage());
        }
    }
    
    // Funções auxiliares para os testes
    public static double[][] multiplicarMatrizes(double[][] A, double[][] B) {
        int n = A.length;
        double[][] C = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double soma = 0;
                for (int k = 0; k < n; k++) {
                    soma += A[i][k] * B[k][j];
                }
                C[i][j] = soma;
            }
        }
        return C;
    }
    
    public static double[][] aplicarPermutacaoLinhas(double[][] A, int[] perm) {
        int n = A.length;
        double[][] P_A = new double[n][n];
        for (int i = 0; i < n; i++) {
            P_A[i] = A[perm[i]].clone();
        }
        return P_A;
    }
    
    public static double[][] aplicarPermutacaoLinhasColunas(double[][] A, int[] permLinhas, int[] permColunas) {
        int n = A.length;
        double[][] PAQ = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                PAQ[i][j] = A[permLinhas[i]][permColunas[j]];
            }
        }
        return PAQ;
    }
}