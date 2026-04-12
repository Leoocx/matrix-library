/**
 * Classe que realiza a decomposição LU (sem pivotamento) de uma matriz A
 * e aplica as mesmas operações elementares a um vetor b.
 * 
 * A decomposição é: A = L * U, onde:
 * - L é triangular inferior com diagonal unitária
 * - U é triangular superior
 * 
 * O vetor b é transformado para b_modificado, que corresponde ao vetor
 * b após as mesmas eliminações.
 * 
 * Esta classe NÃO resolve o sistema; apenas fornece L, U e b_modificado.
 */
public class LUSemPivotamento {
    private double[][] L;
    private double[][] U;
    private double[] bMod;
    private int n;

    public LUSemPivotamento(double[][] A, double[] b) {
        n = A.length;
        // Verificações
        for (int i = 0; i < n; i++)
            if (A[i].length != n)
                throw new IllegalArgumentException("Matriz não quadrada");
        if (b.length != n)
            throw new IllegalArgumentException("b tamanho incorreto");

        // Copia A para U (será modificada)
        U = new double[n][n];
        for (int i = 0; i < n; i++)
            System.arraycopy(A[i], 0, U[i], 0, n);

        bMod = b.clone();
        L = new double[n][n];
        for (int i = 0; i < n; i++)
            L[i][i] = 1.0; // diagonal unitária

        // Eliminação de Gauss (sem pivotamento)
        for (int k = 0; k < n - 1; k++) {
            if (U[k][k] == 0)
                throw new RuntimeException("Pivô zero na coluna " + k + " - matriz singular ou precisa pivotamento");
            for (int i = k + 1; i < n; i++) {
                double fator = U[i][k] / U[k][k];
                L[i][k] = fator;
                for (int j = k; j < n; j++) {
                    U[i][j] -= fator * U[k][j];
                }
                bMod[i] -= fator * bMod[k];
            }
        }
        if (U[n-1][n-1] == 0)
            throw new RuntimeException("Matriz singular");
    }

    public double[][] getL() { return copiarMatriz(L); }
    public double[][] getU() { return copiarMatriz(U); }
    public double[] getBModificado() { return bMod.clone(); }

    private double[][] copiarMatriz(double[][] M) {
        double[][] copia = new double[n][n];
        for (int i = 0; i < n; i++)
            System.arraycopy(M[i], 0, copia[i], 0, n);
        return copia;
    }

    public static void main(String[] args) {
        double[][] A = {
            {1, 2, 1, 0},
            {2, 5, 1, 1},
            {1, 1, 3, 2},
            {0, 1, 2, 6}
        };
        double[] x = {1,1,1,1};

        double[] b = multiplicarMatrizVetor(A, x);


        System.out.println("b:");
            for (double val : b) {
                System.out.printf("%.2f ", val);
            }System.out.println();

        LUSemPivotamento lu = new LUSemPivotamento(A, b);
            System.out.println("L (triangular inferior):");
        imprimirMatriz(lu.getL());
            System.out.println("\nU (triangular superior):");
        imprimirMatriz(lu.getU());
            System.out.println("\nb modificado:");
        imprimirVetor(lu.getBModificado());
    
        /*double[][] LU =multiplicarMatrizes(lu.getL(), lu.getU()) ;
        System.out.println("Matriz LU:");
        imprimirMatriz(LU);*/
    }

    private static void imprimirMatriz(double[][] M) {
        for (double[] linha : M) {
            for (double v : linha) System.out.printf("%8.4f ", v);
            System.out.println();
        }
    }
    private static void imprimirVetor(double[] v) {
        for (double val : v) System.out.printf("%8.4f ", val);
        System.out.println();
    }

    /**
 * Multiplica uma matriz A (n x m) por um vetor x (m) resultando em b (n).
 * 
 * @param A matriz de coeficientes (n linhas, m colunas)
 * @param x vetor de entrada (tamanho m)
 * @return vetor b = A * x (tamanho n)
 * @throws IllegalArgumentException se as dimensões forem incompatíveis
 */
public static double[] multiplicarMatrizVetor(double[][] A, double[] x) {
    // Validações
    if (A == null || x == null)
        throw new IllegalArgumentException("Matriz e vetor não podem ser nulos.");
    
    int n = A.length;          // número de linhas
    if (n == 0)
        throw new IllegalArgumentException("Matriz não pode ter zero linhas.");
    
    int m = A[0].length;       // número de colunas (assumindo matriz retangular)
    // Verifica se todas as linhas têm o mesmo número de colunas
    for (int i = 1; i < n; i++) {
        if (A[i].length != m)
            throw new IllegalArgumentException("Matriz irregular: linha " + i + " tem " + A[i].length + " colunas, esperado " + m);
    }
    
    if (x.length != m)
        throw new IllegalArgumentException("Tamanho do vetor x (" + x.length + ") não coincide com número de colunas da matriz (" + m + ")");
    
    double[] b = new double[n];
    
    // Multiplicação: b[i] = soma_{j=0}^{m-1} A[i][j] * x[j]
    for (int i = 0; i < n; i++) {
        double soma = 0.0;
        for (int j = 0; j < m; j++) {
            soma += A[i][j] * x[j];
        }
        b[i] = soma;
    }
    return b;
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
}