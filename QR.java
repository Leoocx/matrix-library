/**
 * Métodos para ortogonalização de bases (Gram-Schmidt) e decomposição QR de matrizes.
 * 
 */
public class QR {

    /**
     * Ortogonaliza uma base de vetores usando o processo de Gram-Schmidt clássico.
     * 
     * @param base vetores de entrada (linearmente independentes) organizados em uma matriz
     *             onde cada coluna é um vetor da base (dimensão n x m, n = dimensão do espaço,
     *             m = número de vetores, com m ≤ n)
     * @return uma matriz com os vetores ortogonais (colunas), ainda não normalizados
     */
    public static double[][] gramSchmidtOrtogonal(double[][] base) {
        int n = base.length;      // dimensão do espaço
        int m = base[0].length;   // número de vetores na base
        
        double[][] u = new double[n][m]; // vetores ortogonais
        
        for (int j = 0; j < m; j++) {
            // Inicializa u_j como o vetor original
            for (int i = 0; i < n; i++) {
                u[i][j] = base[i][j];
            }
            
            // Subtrai as projeções sobre os vetores ortogonais já calculados (u_0 ... u_{j-1})
            for (int k = 0; k < j; k++) {
                // Calcula o produto escalar <base_j, u_k>
                double prodEscalar = 0.0;
                for (int i = 0; i < n; i++) {
                    prodEscalar += base[i][j] * u[i][k];
                }
                // Calcula a norma ao quadrado de u_k
                double norma2 = 0.0;
                for (int i = 0; i < n; i++) {
                    norma2 += u[i][k] * u[i][k];
                }
                double fator = prodEscalar / norma2;
                // Subtrai a projeção: u_j = u_j - fator * u_k
                for (int i = 0; i < n; i++) {
                    u[i][j] -= fator * u[i][k];
                }
            }
        }
        return u;
    }

    /**
     * Ortonormaliza uma base (vetores ortogonais e unitários).
     * 
     * @param base vetores de entrada (linearmente independentes)
     * @return matriz com vetores ortonormais (colunas)
     */
    public static double[][] gramSchmidtOrtonormal(double[][] base) {
        double[][] u = gramSchmidtOrtogonal(base);
        int n = u.length;
        int m = u[0].length;
        double[][] q = new double[n][m];
        
        // Normaliza cada vetor
        for (int j = 0; j < m; j++) {
            double norma = 0.0;
            for (int i = 0; i < n; i++) {
                norma += u[i][j] * u[i][j];
            }
            norma = Math.sqrt(norma);
            if (norma == 0) {
                throw new RuntimeException("Vetor nulo encontrado (base não é linearmente independente).");
            }
            for (int i = 0; i < n; i++) {
                q[i][j] = u[i][j] / norma;
            }
        }
        return q;
    }

    /**
     * Decomposição QR de uma matriz A (n x m, com n ≥ m e posto m) usando Gram-Schmidt.
     * Retorna Q (n x m) ortogonal (colunas ortonormais) e R (m x m) triangular superior
     * tais que A = Q * R.
     * 
     * @param A matriz de entrada (n x m)
     * @return um array contendo [Q, R]
     */
    public static double[][][] decomposicaoQR(double[][] A) {
        int n = A.length;
        int m = A[0].length;
        if (n < m) {
            throw new IllegalArgumentException("A matriz deve ter número de linhas ≥ número de colunas.");
        }
        
        double[][] Q = new double[n][m];
        double[][] R = new double[m][m];
        
        // Gram-Schmidt modificado (mais estável numericamente)
        // Inicializa Q com cópia de A
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                Q[i][j] = A[i][j];
            }
        }
        
        for (int j = 0; j < m; j++) {
            // Calcula R[j][j] = norma da coluna j de Q (antes da normalização)
            double norma = 0.0;
            for (int i = 0; i < n; i++) {
                norma += Q[i][j] * Q[i][j];
            }
            norma = Math.sqrt(norma);
            if (norma == 0) {
                throw new RuntimeException("Matriz A é singular (coluna " + j + " é combinação linear das anteriores).");
            }
            R[j][j] = norma;
            
            // Normaliza a coluna j de Q
            for (int i = 0; i < n; i++) {
                Q[i][j] /= norma;
            }
            
            // Para as colunas seguintes (k = j+1 .. m-1), subtrai a projeção sobre a coluna j
            for (int k = j + 1; k < m; k++) {
                // Calcula R[j][k] = Q[:,j] · Q[:,k] (vetor ainda não processado)
                double prodEscalar = 0.0;
                for (int i = 0; i < n; i++) {
                    prodEscalar += Q[i][j] * Q[i][k];
                }
                R[j][k] = prodEscalar;
                // Subtrai a projeção: Q[:,k] = Q[:,k] - R[j][k] * Q[:,j]
                for (int i = 0; i < n; i++) {
                    Q[i][k] -= R[j][k] * Q[i][j];
                }
            }
        }
        
        return new double[][][]{Q, R};
    }

    // --------------------------------------------------------------
    // Exemplo de uso: base em R^5 e decomposição QR
    // --------------------------------------------------------------
    public static void main(String[] args) {
        // 1. Criar uma base linearmente independente em R^5 com 5 vetores
        double[][] base = {
            {1, 0, 1, 0, 1},   // vetor 1 (coluna 0)
            {1, 1, 0, 0, 2},   // vetor 2 (coluna 1)
            {0, 1, 1, 1, 0},   // vetor 3 (coluna 2)
            {1, 0, 0, 1, 1},   // vetor 4 (coluna 3)
            {0, 0, 1, 1, 1}    // vetor 5 (coluna 4)
        };
        // A matriz base está no formato: linhas = dimensão (5), colunas = vetores (5)
        // Transpor para que cada coluna seja um vetor (mais comum em álgebra linear)
        double[][] vetores = new double[5][5];
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                vetores[i][j] = base[i][j];
            }
        }
        
        System.out.println("Base original (colunas):");
        imprimirMatriz(vetores);
        
        // Ortogonalização (Gram-Schmidt)
        double[][] ortogonal = gramSchmidtOrtogonal(vetores);
        System.out.println("\nBase ortogonal (não normalizada):");
        imprimirMatriz(ortogonal);
        
        // Ortonormalização
        double[][] ortonormal = gramSchmidtOrtonormal(vetores);
        System.out.println("\nBase ortonormal:");
        imprimirMatriz(ortonormal);
        
        // Verificar ortogonalidade: produto escalar entre colunas i e j deve ser 0 (para i≠j) e 1 (i=j)
        System.out.println("\nVerificação da ortonormalidade (Q^T * Q):");
        double[][] QTQ = multiplicarMatrizes(transpor(ortonormal), ortonormal);
        imprimirMatriz(QTQ);
        
        // 2. Decomposição QR de uma matriz não singular
        // Exemplo: matriz 5x5 com determinante não nulo
        double[][] A = {
            {2, 1, 0, 1, 3},
            {1, 2, 1, 0, 1},
            {0, 1, 3, 2, 0},
            {1, 0, 2, 4, 1},
            {3, 1, 0, 1, 2}
        };
        
        System.out.println("\n\nMatriz A (não singular):");
        imprimirMatriz(A);
        
        double[][][] QR = decomposicaoQR(A);
        double[][] Q = QR[0];
        double[][] R = QR[1];
        
        System.out.println("\nMatriz Q (ortogonal):");
        imprimirMatriz(Q);
        System.out.println("\nMatriz R (triangular superior):");
        imprimirMatriz(R);
        
        // Verificação: A = Q * R
        double[][] QRproduto = multiplicarMatrizes(Q, R);
        System.out.println("\nVerificação A = Q * R:");
        imprimirMatriz(QRproduto);
        
        // Erro máximo
        double maxErro = 0.0;
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[0].length; j++) {
                maxErro = Math.max(maxErro, Math.abs(A[i][j] - QRproduto[i][j]));
            }
        }
        System.out.printf("\nErro máximo da decomposição: %.2e\n", maxErro);
    }
    
    // Funções auxiliares para manipulação de matrizes
    private static double[][] transpor(double[][] M) {
        int linhas = M.length;
        int colunas = M[0].length;
        double[][] Mt = new double[colunas][linhas];
        for (int i = 0; i < linhas; i++) {
            for (int j = 0; j < colunas; j++) {
                Mt[j][i] = M[i][j];
            }
        }
        return Mt;
    }
    
    private static double[][] multiplicarMatrizes(double[][] A, double[][] B) {
        int n = A.length;
        int m = A[0].length;
        int p = B[0].length;
        double[][] C = new double[n][p];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < p; j++) {
                double soma = 0.0;
                for (int k = 0; k < m; k++) {
                    soma += A[i][k] * B[k][j];
                }
                C[i][j] = soma;
            }
        }
        return C;
    }
    
    private static void imprimirMatriz(double[][] M) {
        for (double[] linha : M) {
            for (double val : linha) {
                System.out.printf("%8.4f ", val);
            }
            System.out.println();
        }
    }
}