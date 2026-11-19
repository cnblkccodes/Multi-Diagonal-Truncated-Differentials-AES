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

double compute_triple_sum2(int j) {
double total2 = 0.0;
    for (int m1 = 1; m1 <= 4; m1++) {
            for (int m2 = m1; m2 <= 4; m2++) {
                double term2 = (2*delta(m2-m1) + delta(m1-m2+1)) * combination(4,m1) * pow(2,31) * pow(pow(2,8)-1,m1) * combination(4,m2) * pow(2,31) * pow(pow(2,8)-1,m2) * compute_triple_sum(m1,m2,j) * pow(2,-62) * pow(pow(2,32)-1,-2);
                total2 += term2;

            }
    }
return total2;
}




double compute_triple_sum3(int i, int j) {
    double total3 = 0.0;
    int min_ji = (4-j < i) ? 4-j : i;
    for (int l = 0; l <= min_ji-1; l++) {
    double term3 = combination(4, l) * P_ij(4-j,l);
        for (int t = 0; t <= i-l-1; t++) {
            int sign2 = (t % 2 == 0) ? 1 : -1;

            double term2 = term3 * combination(i,t) * sign2 * pow(P_igj(4 - l, 4-i+t), 4);
            total3 += term2;
        }
    }
    return total3;
}

double compute_triple_sum4(int i,int j) {
    double total4 = 0.0;
    total4 = compute_triple_sum3(i,j) * combination(4,j) * pow(1-pow(2,-8),4-j) * pow(2,-32*i-8*j+32) * pow(1-pow(2,-32),-i);
    return total4;
}


int main() {
    int m1,m2,i,j;
    printf("Enter i: ");
    scanf("%d", &i);
    printf("Enter j: ");
    scanf("%d", &j);

    double triple_sum_result1 = combination(4,j) * compute_triple_sum2(j);
    double triple_sum_result = compute_triple_sum(m1,m2,j);
    double triple_sum_result2 = 4 * compute_triple_sum4(i,j);


    printf("\n Average * combination(4,m) = Pr[2->j] * 4");
    printf("\n %.20f =  %.20f", log2(triple_sum_result1), log2(triple_sum_result2));
    return 0;
}

