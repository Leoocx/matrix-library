class SubstituicaoDireta{

    /**
 * Resolve um sistema Lx = b, onde L é uma matriz triangular inferior.
 * 
 * Algoritmo de substituição direta: começa da primeira equação (linha 0) e
 * desce até a última, calculando x[i] a partir das variáveis já conhecidas x[0...i-1].
 * 
 * @param L Matriz triangular inferior (n x n). A diagonal não pode conter zeros.
 * @param b Vetor dos termos independentes (tamanho n)
 * @return Vetor solução x (tamanho n)
 * @throws IllegalArgumentException se as dimensões são inválidas ou matriz é singular
 */
public static double[] substituicaoDireta(double[][] L, double[] b) {
    int n = L.length;
    
    // Verificação de consistência dimensional
    if (b.length != n) {
        throw new IllegalArgumentException("Dimensões incompatíveis: L é " + n + "x" + n +
" mas b tem " + b.length + " elementos");
    }
    
    double[] x = new double[n];
    
    // Itera de cima para baixo (da primeira linha para a última)
    for (int i = 0; i < n; i++) {
        double soma = 0.0;
        
        // Acumula a contribuição das variáveis já calculadas (x[j] com j < i)
        for (int j = 0; j < i; j++) {
            soma += L[i][j] * x[j];
        }
        
        // Calcula x[i] = (b[i] - soma) / L[i][i]
        double pivô = L[i][i];
        if (pivô == 0) {
            throw new IllegalArgumentException("Matriz singular (elemento diagonal zero na linha " + i + ")");
        }
        x[i] = (b[i] - soma) / pivô;
    }
    
    return x;
}

    public static void main(String[] args) {
        double[][] L = {
            {1,  0,  0, 0},
            {0,  1,  0, 0},
            {4, 3,  1, 0},
            {3, 4, 1, 1}
        };
        double[] b2 = {3, 7, 19, 17};
        
        double[] x2 = substituicaoDireta(L, b2);
        System.out.println("\nSolução via substituição direta:");
        System.out.printf("x = %.2f, y = %.2f, w = %.2f, z = %.2f\n", x2[0], x2[1], x2[2], x2[3]);
        
    }
    
    


}