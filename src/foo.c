#include <stdio.h>

void foo(int *nin, double *x) {
  int n = nin[0];

  for (int i = 0; i < n; ++i) {
    x[i] = x[i] * x[i];
  }
}
