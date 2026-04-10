/**
 * Classe com operações de projeção e ortogonalidade entre vetores.
 * 
 */
public class AlgebraVetorial {

    /**
     * Calcula o produto escalar (dot product) entre dois vetores.
     * @param u primeiro vetor
     * @param v segundo vetor
     * @return u · v
     */
    public static double produtoEscalar(double[] u, double[] v) {
        if (u.length != v.length) {
            throw new IllegalArgumentException("Vetores devem ter mesma dimensão.");
        }
        double soma = 0.0;
        for (int i = 0; i < u.length; i++) {
            soma += u[i] * v[i];
        }
        return soma;
    }

    /**
     * Calcula a norma (magnitude) de um vetor.
     * @param v vetor
     * @return ||v||
     */
    public static double norma(double[] v) {
        double soma = 0.0;
        for (double comp : v) {
            soma += comp * comp;
        }
        return Math.sqrt(soma);
    }

    /**
     * Projeção escalar de u sobre v (também chamada de componente escalar).
     * Resultado: comp_v(u) = (u·v) / ||v||
     * É o comprimento (com sinal) da projeção vetorial.
     * 
     * @param u vetor que será projetado
     * @param v vetor sobre o qual se projeta (não pode ser nulo)
     * @return componente escalar (double)
     */
    public static double projecaoEscalar(double[] u, double[] v) {
        double normaV = norma(v);
        if (normaV == 0) {
            throw new IllegalArgumentException("Vetor v não pode ser nulo.");
        }
        return produtoEscalar(u, v) / normaV;
    }

    /**
     * Projeção vetorial de u sobre v.
     * Resultado: proj_v(u) = ((u·v)/(v·v)) * v
     * É o vetor que representa a sombra de u na direção de v.
     * 
     * @param u vetor que será projetado
     * @param v vetor sobre o qual se projeta (não pode ser nulo)
     * @return vetor projeção (mesma dimensão de u e v)
     */
    public static double[] projecaoVetorial(double[] u, double[] v) {
        double produtoVV = produtoEscalar(v, v);
        if (produtoVV == 0) {
            throw new IllegalArgumentException("Vetor v não pode ser nulo.");
        }
        double fator = produtoEscalar(u, v) / produtoVV;
        double[] proj = new double[v.length];
        for (int i = 0; i < v.length; i++) {
            proj[i] = fator * v[i];
        }
        return proj;
    }

    /**
     * Retorna um vetor não nulo que seja ortogonal ao vetor dado.
     * 
     * Este é um exemplo de como obter um vetor do complemento ortogonal
     * (espaço de todos os vetores w tais que v·w = 0).
     * O complemento ortogonal de um vetor não nulo em R^n tem dimensão n-1.
     * 
     * Algoritmo: Encontra a primeira componente não nula de v e constrói w
     * que tenha um 1 na posição seguinte (circularmente) e zeros nas demais,
     * ajustado para garantir ortogonalidade.
     * 
     * @param v vetor de entrada (não pode ser nulo)
     * @return vetor w tal que v·w = 0 e w não nulo
     */
    public static double[] complementoOrtogonal(double[] v) {
        // Verifica se v é nulo
        boolean nulo = true;
        for (double comp : v) {
            if (comp != 0.0) {
                nulo = false;
                break;
            }
        }
        if (nulo) {
            throw new IllegalArgumentException("Vetor nulo não possui complemento ortogonal único.");
        }

        int n = v.length;
        double[] w = new double[n];
        
        // Encontra o índice da primeira componente não nula
        int idx = 0;
        while (idx < n && v[idx] == 0) {
            idx++;
        }
        // Agora v[idx] != 0
        // Estratégia: w será zero em todas posições exceto em idx e (idx+1 mod n)
        int next = (idx + 1) % n;
        w[next] = 1.0;
        // Para garantir ortogonalidade: v[idx]*w[idx] + v[next]*w[next] = 0
        // Como w[next]=1, temos v[idx]*w[idx] + v[next] = 0 => w[idx] = -v[next]/v[idx]
        w[idx] = -v[next] / v[idx];
        // As demais componentes permanecem 0
        return w;
    }

    public static void main(String[] args) {
        // Exemplo 1: Vetores no R²
        double[] u = {3, 4};
        double[] v = {1, 0};
        
        System.out.println("Exemplo em R²:");
        System.out.printf("u = (%.1f, %.1f), v = (%.1f, %.1f)\n", u[0], u[1], v[0], v[1]);
        
        double esc = projecaoEscalar(u, v);
        System.out.printf("Projeção escalar de u sobre v: %.2f\n", esc);
        
        double[] vet = projecaoVetorial(u, v);
        System.out.printf("Projeção vetorial de u sobre v: (%.2f, %.2f)\n", vet[0], vet[1]);
        
        double[] w = complementoOrtogonal(v);
        System.out.printf("Vetor ortogonal a v: (%.2f, %.2f)\n", w[0], w[1]);
        System.out.printf("Verificação v·w = %.2f\n\n", produtoEscalar(v, w));
        
        // Exemplo 2: Vetores no R³
        double[] a = {2, -1, 3};
        double[] b = {4, 2, 0};
        
        System.out.println("Exemplo em R³:");
        System.out.printf("a = (%.1f, %.1f, %.1f), b = (%.1f, %.1f, %.1f)\n", 
a[0], a[1], a[2], b[0], b[1], b[2]);
        
        double[] proj = projecaoVetorial(a, b);
        System.out.printf("Projeção vetorial de a sobre b: (%.2f, %.2f, %.2f)\n", 
proj[0], proj[1], proj[2]);
        
        double[] ort = complementoOrtogonal(b);
        System.out.printf("Vetor ortogonal a b: (%.2f, %.2f, %.2f)\n", ort[0], ort[1], ort[2]);
        System.out.printf("Verificação b·ort = %.2f\n", produtoEscalar(b, ort));
    }
}