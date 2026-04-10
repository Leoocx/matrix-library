/**
 * Resolução de sistemas lineares utilizando o método de eliminação de Gauss-Jordan.
 * 
 * O método transforma a matriz aumentada [A | b] em uma matriz onde a parte de A
 * se torna a matriz identidade, resultando diretamente no vetor solução x no lugar de b.
 * 
 * Complexidade: O(n³) para uma matriz n x n.
 * 
 * Vantagem: Não requer substituição retroativa (já que a matriz se torna identidade).
 * Desvantagem: Cerca de 50% mais operações que a eliminação de Gauss simples.
 */
public class GaussJordan {

    /**
     * Resolve o sistema linear A * x = b usando o método de Gauss-Jordan.
     * 
     * @param A Matriz dos coeficientes (n x n). Será modificada.
     * @param b Vetor dos termos independentes (tamanho n). Será modificado.
     * @return Vetor solução x (tamanho n), ou null se a matriz for singular.
     */
    public static double[] resolver(double[][] A, double[] b) {
        int n = A.length;
        
        // Cria uma cópia da matriz aumentada [A | b] para não modificar os originais?
        // (opcional: se quiser preservar A e b, descomente as linhas abaixo)
        // double[][] M = new double[n][n+1];
        // for (int i = 0; i < n; i++) {
        //     System.arraycopy(A[i], 0, M[i], 0, n);
        //     M[i][n] = b[i];
        // }
        // Mas aqui vamos trabalhar diretamente com A e b por eficiência.
        
        // Etapa 1: Escalonamento para obter matriz identidade
        for (int k = 0; k < n; k++) {
            // --- Pivotamento parcial: encontrar o maior pivô (em módulo) na coluna k
            int maxIndex = k;
            double maxVal = Math.abs(A[k][k]);
            for (int i = k + 1; i < n; i++) {
                if (Math.abs(A[i][k]) > maxVal) {
                    maxVal = Math.abs(A[i][k]);
                    maxIndex = i;
                }
            }
            
            // Se o maior pivô é zero, a matriz é singular
            if (maxVal == 0) {
                System.err.println("Matriz singular! Sistema não tem solução única.");
                return null;
            }
            
            // Trocar linha k pela linha maxIndex, se necessário
            if (maxIndex != k) {
                double[] tempRow = A[k];
                A[k] = A[maxIndex];
                A[maxIndex] = tempRow;
                
                double tempB = b[k];
                b[k] = b[maxIndex];
                b[maxIndex] = tempB;
            }
            
            // --- Normalizar a linha do pivô (dividir toda a linha pelo pivô)
            double pivô = A[k][k];
            for (int j = k; j < n; j++) {
                A[k][j] /= pivô;
            }
            b[k] /= pivô;
            
            // --- Eliminar as outras linhas (tornar zero os elementos da coluna k em todas as outras linhas)
            for (int i = 0; i < n; i++) {
                if (i != k) {
                    double fator = A[i][k];
                    if (fator != 0) {
                        // Atualizar a linha i: Li <- Li - fator * Lk
                        for (int j = k; j < n; j++) {
                            A[i][j] -= fator * A[k][j];
                        }
                        b[i] -= fator * b[k];
                    }
                }
            }
        }
        
        // Após o processo, a matriz A se tornou identidade e b contém a solução.
        return b;
    }
    
    public static void main(String[] args) {
        // Exemplo 1: Sistema 3x3
        double[][] A1 = {
            {2, 1, -1},
            {1, 2, 1},
            {1, 1, 1}
        };
        double[] b1 = {-3, 3, 2};
        
        double[] x1 = resolver(A1, b1);
        if (x1 != null) {
            System.out.println("Solução do sistema 3x3:");
            for (int i = 0; i < x1.length; i++) {
                System.out.printf("x%d = %.6f\n", i+1, x1[i]);
            }
        }
        
        // Exemplo 2: Sistema 2x2
        double[][] A2 = {{2, 3}, {4, 9}};
        double[] b2 = {6, 15};
        double[] x2 = resolver(A2, b2);
        if (x2 != null) {
            System.out.println("\nSolução do sistema 2x2:");
            for (int i = 0; i < x2.length; i++) {
                System.out.printf("x%d = %.6f\n", i+1, x2[i]);
            }
        }
        
        // Exemplo 3: Sistema 4x4
        double[][] A3 = {
            {4, 1, 2, -3},
            {-3, 3, -1, 4},
            {-1, 2, 5, 1},
            {5, 4, 3, -1}
        };
        double[] b3 = {-16, 20, -4, -10};
        double[] x3 = resolver(A3, b3);
        if (x3 != null) {
            System.out.println("\nSolução do sistema 4x4:");
            for (int i = 0; i < x3.length; i++) {
                System.out.printf("x%d = %.6f\n", i+1, x3[i]);
            }
        }
    }
}