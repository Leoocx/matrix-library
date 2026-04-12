/**
 * Classe que implementa o método de eliminação de Gauss para resolução de sistemas lineares.
 * A solução é obtida para sistemas da forma A * x = b, onde A é uma matriz quadrada de coeficientes
 * e b é o vetor de termos independentes.
 * 
 * O algoritmo utiliza pivotamento parcial (busca pelo maior pivô em módulo) para melhorar a estabilidade numérica.
 * Complexidade de tempo: O(n³), onde n é a ordem do sistema.
 */
public class GaussSolver {    

    /**
     * Resolve um sistema linear utilizando eliminação de Gauss com pivotamento parcial.
     * 
     * Etapas do método:
     * 1. Escalonamento (triangularização) da matriz aumentada [A | b].
     * 2. Substituição retroativa para obter a solução x.
     * 
     * @param A Matriz quadrada dos coeficientes (n x n). Será modificada durante o escalonamento.
     * @param b Vetor dos termos independentes (tamanho n). Também será modificado.
     * @return Vetor com a solução do sistema (tamanho n), ou null se a matriz for singular (pivô zero).
     */
    public static double[] gaussSolver(double[][] A, double[] b) {
        
        // ------------------------------------------------------------
        // ETAPA 1: ESCALONAMENTO (triangularização superior)
        // Itera sobre cada coluna k (pivô) até a penúltima linha.
        // ------------------------------------------------------------
        for (int k = 0; k < A.length - 1; k++) {
            
            // --- Pivotamento parcial: encontrar o maior elemento (em módulo) na coluna k,
            // a partir da linha k (inclusive).
            double max = Math.abs(A[k][k]);   // Assume o pivô atual como o maior
            int maxIndex = k;                 // Índice da linha que contém o maior elemento
            
            // Percorre as linhas abaixo de k para buscar um elemento de maior módulo
            for (int i = k + 1; i < A.length; i++) {
                if (max < Math.abs(A[i][k])) {
                    max = Math.abs(A[i][k]);
                    maxIndex = i;            // Nova linha com maior pivô
                }
            }
            
            // Se a linha com o maior elemento não for a linha atual (k), troca as linhas inteiras
            // tanto na matriz A quanto no vetor b.
            if (maxIndex != k) {
                // Troca as linhas k e maxIndex na matriz A
                for (int j = 0; j < A.length; j++) {
                    double temp = A[k][j];
                    A[k][j] = A[maxIndex][j];
                    A[maxIndex][j] = temp;
                }
                // Troca os respectivos elementos no vetor b
                double temp = b[k];
                b[k] = b[maxIndex];
                b[maxIndex] = temp;
            }
            
            // Após o pivotamento, verifica se o pivô (agora em A[k][k]) é zero.
            // Se for zero, a matriz é singular (determinante zero) e o sistema não tem solução única.
            if (A[k][k] == 0) {
                return null;   // Sistema impossível ou indeterminado
            } 
            else {
                // --- Eliminação: zerar os elementos abaixo do pivô na coluna k
                // Percorre todas as linhas abaixo da linha do pivô
                for (int m = k + 1; m < A.length; m++) {
                    // Calcula o fator que irá zerar A[m][k] usando o pivô A[k][k]
                    // F = - (elemento a ser zerado) / (pivô)
                    double F = -A[m][k] / A[k][k];
                    
                    // Como A[m][k] será zerado, já o colocamos como 0 explicitamente (evita uma iteração no loop abaixo)
                    A[m][k] = 0;
                    
                    // Atualiza o termo independente b[m] de acordo com a operação de linha: L_m <- L_m + F * L_k
                    b[m] = b[m] + F * b[k];
                    
                    // Atualiza os demais elementos da linha m (colunas à direita da coluna k)
                    // Começa em j = k+1 porque os elementos à esquerda (colunas 0..k) já estão zerados ou não precisam ser alterados.
                    for (int l = k + 1; l < A.length; l++) {
                        A[m][l] = A[m][l] + F * A[k][l];
                    }
                }
            }
        }
        
        // ------------------------------------------------------------
        // ETAPA 2: SUBSTITUIÇÃO RETROATIVA
        // Após o escalonamento, a matriz A é triangular superior.
        // Resolve-se de baixo para cima: x[n-1], x[n-2], ..., x[0].
        // ------------------------------------------------------------
        double[] X = new double[A.length];   // Vetor que armazenará a solução
        
        // Itera das últimas linhas para as primeiras
        for (int i = A.length - 1; i >= 0; i--) {
            // Inicializa X[i] com o valor de b[i] (termo independente após as transformações)
            X[i] = b[i];
            
            // Subtrai as contribuições das variáveis já calculadas (X[j] para j > i)
            for (int j = i + 1; j < A.length; j++) {
                X[i] = X[i] - X[j] * A[i][j];
            }
            
            // Divide pelo coeficiente diagonal A[i][i] (pivô) para obter o valor da variável i
            X[i] = X[i] / A[i][i];
        }
        
        return X;   // Retorna o vetor solução
    }

    public static void main(String args[]) {
        
        // Teste 1: Sistema 3x3
        double A[][] = {{2, 1, -1}, {1, 2, 1}, {1, 1, 1}};
        double b[] = {-3, 3, 2};
        double[] x = gaussSolver(A, b);
        System.out.printf("x1 = %f\nx2 = %f\nx3 = %f\n", x[0], x[1], x[2]);

        // Teste 2: Outro sistema 3x3
        double A2[][] = {{1, 1, 1}, {-2, 1, 1}, {1, 3, 1}};
        double b2[] = {2, 5, 4};
        x = gaussSolver(A2, b2);
        System.out.printf("\nx1 = %f\nx2 = %f\nx3 = %f\n", x[0], x[1], x[2]);
        
        // Teste 3: Sistema 3x3 com coeficientes decimais
        double A3[][] = {{3, 2, -1}, {2, -2, 4}, {-1, 0.5, -1}};
        double b3[] = {1, -2, 0};
        x = gaussSolver(A3, b3);
        System.out.printf("\nx1 = %f\nx2 = %f\nx3 = %f\n", x[0], x[1], x[2]);
        
        // Teste 4: Sistema 2x2
        double A4[][] = {{2,3}, {4,9}};
        double b4[] = {6,15};
        x = gaussSolver(A4, b4);
        System.out.printf("\nx1 = %f\nx2 = %f\n", x[0], x[1]);

        // Teste 5: Sistema 4x4
        double A5[][] = 
        {{2, 1, 1, 0}, 
        {0, 1, 1, 1},
        {8, 7, 9, 5}, 
        {6, 7, 9, 8}};
        double b5[] = {3, 7, 19, 17};
        x = gaussSolver(A5, b5);
        System.out.printf("\nx1 = %f\nx2 = %f\nx3 = %f\nx4 = %f\n", x[0], x[1], x[2], x[3]);

        // Teste 6: Sistema 5x5
        double A6[][] = {{4, 1, 2, -3, 5}, {-3, 3, -1, 4, -2}, {-1, 2, 5, 1, 3}, {5, 4, 3, -1, 2}, {1, -2, 3, -4, 5}};
        double b6[] = {-16, 20, -4, -10, 3};
        x = gaussSolver(A6, b6);
        System.out.printf("\nx1 = %f\nx2 = %f\nx3 = %f\nx4 = %f\nx5 = %f\n", x[0], x[1], x[2], x[3], x[4]);
    }
}