class RetroSubstituicao{

    /**
 * Resolve um sistema linear Ux = b, onde U é uma matriz triangular superior.
 * 
 * Algoritmo de retrosubstituição: começa da última equação (linha n-1) e
 * sobe até a primeira, calculando x[i] a partir das variáveis já conhecidas x[i+1...n-1].
 * 
 * @param U Matriz triangular superior (n x n). A diagonal não pode conter zeros.
 * @param b Vetor dos termos independentes (tamanho n)
 * @return Vetor solução x (tamanho n)
 * @throws IllegalArgumentException se as dimensões são inválidas ou matriz é singular
 */
public static double[] retrosubstituicao(double[][] U, double[] b) {
    int n = U.length;
    
    // Verificação de consistência dimensional
    if (b.length != n) {
        throw new IllegalArgumentException("Dimensões incompatíveis: U é " + n + "x" + n +
" mas b tem " + b.length + " elementos");
    }
    
    double[] x = new double[n];
    
    // Itera de baixo para cima (da última linha para a primeira)
    for (int i = n - 1; i >= 0; i--) {
        double soma = 0.0;
        
        // Acumula a contribuição das variáveis já calculadas (x[j] com j > i)
        for (int j = i + 1; j < n; j++) {
            soma += U[i][j] * x[j];
        }
        
        // Calcula x[i] = (b[i] - soma) / U[i][i]
        double pivô = U[i][i];
        if (pivô == 0) {
            throw new IllegalArgumentException("Matriz singular (elemento diagonal zero na linha " + i + ")");
        }
        x[i] = (b[i] - soma) / pivô;
    }
    
    return x;
}

/*


    Sistema triangular superior
    2x + 3y -  z = 8
    4y + 2z = 6
    5z = 10
    double[][] U = {
        {2,  3, -1},
        {0,  4,  2},
        {0,  0,  5}
    };
    double[] b1 = {8, 6, 10};
    
    double[] x1 = retrosubstituicao(U, b1);
    System.out.println("Solução via retrosubstituição:");
    System.out.printf("x = %.2f, y = %.2f, z = %.2f\n", x1[0], x1[1], x1[2]);
*/




}