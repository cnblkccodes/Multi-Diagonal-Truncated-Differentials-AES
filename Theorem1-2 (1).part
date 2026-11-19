#include <stdio.h>
#include <math.h>


long double delta(int m) {
    return (m > 0) ? 1.0 : 0.0;
}

int min_(int a, int b) {
    return (a < b) ? a : b;
}

double P_i0(int i) {
    return 1 - 4 * delta(i - 1) * pow(255, -1)
           + 10 * delta(i - 2) * pow(255, -2)
           - 20 * delta(i - 3) * pow(255, -3);
}

double P_ij(int i, int j) {
    if (j > 0) {
        return delta(i - j) * pow(255, -j)
               - 4 * delta(i - j - 1) * pow(255, -(j + 1))
               + 10 * delta(i - j - 2) * pow(255, -(j + 2));
    } else {
        return P_i0(i);
    }
}

double P_igj(int i, int j) {
    return delta(i - j) * pow(255, -j)
           - j * (delta(i - j - 1) * pow(255, -(j + 1))
           - delta(i - j - 2) * pow(255, -(j + 2)));
}

double combination(int n, int k) {
    if (k > n || k < 0) return 0;
    if (k == 0 || k == n) return 1;
    double result = 1;
    for (int i = 1; i <= k; i++) {
        result *= (n - i + 1) / (double)i;
    }
    return result;
}

double permutation(int n, int r) {
    if (r > n || r < 0) return 0;
    double result = 1;
    for (int i = 0; i < r; i++) {
        result *= (n - i);
    }
    return result;
}

double factorial(int n) {
    if (n < 0) return 0;
    if (n == 0 || n == 1) return 1;
    double result = 1;
    for (int i = 2; i <= n; i++) {
        result *= i;
    }
    return result;
}

double compute_modified_double_sum(int m, int i) {
    double total = 0.0;

    int min_mi = (m < 4-i) ? m : 4-i;
    for (int l = 0; l <= min_mi-1; l++) {
        for (int j = 0; j <= 3-l-i; j++) {
            double sign = (j % 2 == 0) ? 1.0 : -1.0;
            double term = combination(4, l) * combination(4, i) * combination(4-i, j) * sign * P_ij(m, l) * pow(P_igj(4 - l, i+j), 4.0);
            total += term;
        }
    }

    return total;
}


double compute_divisor(int i) {
    return combination(4,i) * pow(2, -32*i) * pow(1 - pow(2, -32), 4-i) * pow(1 - pow(2, -128), -1);
}

int main() {
    int m,i;
    printf("Enter m: ");
    scanf("%d", &m);
    printf("Enter i: ");
    scanf("%d", &i);

    double triple_sum_result = compute_modified_double_sum(m,i);


    double divisor = compute_divisor(i);
    double final_result_divided = triple_sum_result / divisor;

    printf("Probability of Theorem1(2.part) for (m,i)=(%d,%d): %.40f\n", m,i,triple_sum_result);
    printf("Probability of Theorem1(2.part) for (m,i)=(%d,%d)(log2): %.25f\n", m,i,log2(triple_sum_result));
    printf("Ratio of Theorem1(2.part) to random permutation for (m,i)=(%d,%d): %.25f\n", m,i,final_result_divided);


    return 0;
}
