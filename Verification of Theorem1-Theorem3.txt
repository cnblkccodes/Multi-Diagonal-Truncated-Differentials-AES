#include <stdio.h>
#include <math.h>

double delta(int m) {
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

double compute_triple_sum1(int i, int m) {
    double total1 = 0.0;
    int min_mi = (4-m < i) ? 4-m : i;
    for (int l = 0; l <= min_mi-1; l++) {
    double term1 = combination(4, l) * P_ij(4-m,l);
        for (int j = 0; j <= i-l-1; j++) {
            int sign = (j % 2 == 0) ? 1 : -1;

            double term2 = term1 * combination(i,j) * sign * pow(P_igj(4 - l, 4-i+j), 4);
            total1 += term2;
        }
    }
    return total1;
}

double compute_triple_sum2(int i,int m) {
    double total = 0.0;
    total = compute_triple_sum1(i,m) * combination(4,m) * pow(1-pow(2,-8),4-m) * pow(2,-32*i-8*m+32) * pow(1-pow(2,-32),-i);
    return total;
}

double compute_modified_double_sum3(int m, int i) {
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


double compute_modified_double_sum4(int j) {
    double total = 0.0;

    for (int m = 1; m <= 4; m++) {
        double term = combination(4, m) * pow(2,31) * pow(pow(2,8)-1,m) * compute_modified_double_sum3(m,3) * combination(4, j) * P_ij(4,j) * pow(4,-1) * pow(2,-31) * pow(pow(2,32)-1,-1);
        total += term;
    }

    return total;
}


int main() {
    int j;
    printf("Enter the number of passive bytes in the active inverse diagonal of the ciphertext pair (denoted by j): ");
    scanf("%d", &j);

    double average= compute_modified_double_sum4(j);
    double resultoftheorem3 = compute_triple_sum2(1,j);


    printf("\n the base-2 logarithm of the probability weighted average with probability of MC_4 = the base-2 logarithm of Pr[1=>j]",j);
    printf("\n %.20f =  %.20f", log2(average), log2(resultoftheorem3));


    return 0;
}
