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

double combination(int n, int k) {
    if (k > n || k < 0) return 0;
    double result = 1;
    for (int i = 1; i <= k; i++) {
        result *= (n - (k - i)) / (double)i;
    }
    return result;
}

double compute_modified_double_sum(int i, int j) {
    double total = 0.0;
        int min_ij = (i < 4-j) ? i : 4-j;
        for (int l = 0; l <= min_ij - 1; l++) {
            double term = 4 * combination(4, l) * pow(P_ij(i, l), 4) *
                          pow(P_ij(4, 3), (4 - l)) *P_ij(4-l, j);
            total += term;
        }
    return total;
}

double compute_divisor(int j) {
    return pow(2, -126+8*(4-j)) * pow(1 - pow(2, -8),4-j) * pow(1 - pow(2, -128), -1);
}

int main() {
    int i, j;

    printf("Enter i: ");
    scanf("%d", &i);

    printf("Enter j: ");
    scanf("%d", &j);

    double double_sum_result = compute_modified_double_sum(i,j);
    double divisor = compute_divisor(j);
    double final_result_divided = double_sum_result / divisor;

    printf("Probability of Theorem2 for (i,j)=(%d,%d): %.40f\n", i,j, double_sum_result);
    printf("Probability of Theorem2 for (i,j)=(%d,%d)(log2): %.25f\n", i,j,log2(double_sum_result));
    printf("Ratio of Theorem2 to random permutation for (i,j)=(%d,%d): %.25f\n", i,j,final_result_divided);


    return 0;
}
