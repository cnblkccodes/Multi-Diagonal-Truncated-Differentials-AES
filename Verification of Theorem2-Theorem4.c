#include <stdio.h>
#include <math.h>

double delta(int m) {
    return (m > 0) ? 1.0 : 0.0;
}

double P_i0(int i) {
    return 1 - 4 * delta(i - 1) * pow(255, -1) + 10 * delta(i - 2) * pow(255, -2) - 20 * delta(i - 3) * pow(255, -3);
}

double P_ij(int i, int j) {
    if (j > 0) {
        return delta(i - j) * pow(255, -j) - 4 * delta(i - j - 1) * pow(255, -(j + 1)) + 10 * delta(i - j - 2) * pow(255, -(j + 2));
    } else {
        return P_i0(i);
    }
}

double P_igj(int i, int j) {
    return delta(i - j) * pow(255, -j) - j * (delta(i - j - 1) * pow(255, -(j + 1)) - delta(i - j - 2) * pow(255, -(j + 2)));
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


double compute_triple_sum(int m1, int m2, int j) {
double total = 0.0;
int min_3m1 = (3 < m1-1) ? 3 : m1-1;
int min_1j = (1 < 5-j) ? 1 : 5-j;
    for (int i1 = 0; i1 <= min_3m1; i1++) {
    int min_4i1m2 = (4-i1 < m2-1) ? 4-i1 : m2-1;
    int min_ss = (min_4i1m2 < 3) ? min_4i1m2 : 3;
            for (int i2 = 0; i2 <= min_ss; i2++) {
                for (int t = 0; t <=min_1j; t++) {
                double term = 4*(delta(1-t)*delta(i1+i2)+delta(1-i1-i2))*combination(4,t)*combination(4,i1)*combination(4-i1,i2)*P_ij(m1,i1)*P_ij(m2,i2)*pow(P_ij(2,t),4-i1-i2)*pow(P_ij(4,3),4-t)*P_ij(4-t,j);
                total += term;
                }
            }
    }
return total;
}

double compute_modified_double_sum2(int i, int j) {
    double total = 0.0;
        int min_ij = (i < 4-j) ? i : 4-j;
        for (int l = 0; l <= min_ij - 1; l++) {
            double term = 4 * combination(4, l) * pow(P_ij(i, l), 4) *
                          pow(P_ij(4, 3), (4 - l)) *P_ij(4-l, j);
            total += term;
        }
    return total;
}

int main() {
    int j;
    printf("Enter the number of passive bytes in the active inverse diagonal of the ciphertext pair (denoted by j): ");
    scanf("%d", &j);


    double resultoftheorem4 = compute_triple_sum(1,1,j);
    double resultoftheorem2 = compute_modified_double_sum2(2,j);


    printf("\n the base-2 logarithm of the probability of Theorem2 for %d = the base-2 logarithm of the probability of Theorem4 for %d",j,j);
    printf("\n %.20f =  %.20f", log2(resultoftheorem2), log2(resultoftheorem4));


    return 0;
}
