#include<iostream>
using namespace std;

#include "ecos.h"

int main () {
  cout<<"\nEnter number of variables/dimensions: ";
  int n;
  cin>>n;
  pfloat *C = new pfloat[n];
  C[0] = 1;
  for(int i=1;i<n;i++) {
    C[i] = C[i-1] + i + 1;
  }

  int m = n+1;
  int p = 0; // No equality constraints.
  int l = m;
  int ncones = 0;
  idxint *q = new idxint[n];
  for (int i=0;i<n;i++) { q[i] = 1; };
  int e = 0;
  // G is size mxn:
  // [ -1 -1 -1 -1 .. -1]
  // [-1 0 0 00 .... 0 0]
  // [0 -1 0 0 0 ... 0 0]
  // .
  // .
  // .
  // [0 0 0 ... 0 0 0 -1]
  // Gpr: -1 -1 ... 2*n times: size 2n
  // Gjc: 0 2 4 6 8 ... 2*n: size n+1
  // Gir: 0 1 0 2 0 .... 0 n: size 2n
  // https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_.28CSC_or_CCS.29
  pfloat* Gpr = new pfloat[2*n];
  idxint* Gjc = new idxint[n+1];
  idxint* Gir = new idxint[2*n];
  for (int i=0;i<2*n;i++) {
    Gpr[i] = -1.0;
    Gir[i] = (i%2==0) ? 0 : (i+1)/2;
    if (i<n+1) {
      Gjc[i] = 2*i;
    }
  }
  // No equality constraints. Hence these A* arrays are NULL.
  pfloat* Apr = NULL;
  idxint* Ajc = NULL;
  idxint* Air = NULL;
  pfloat* h = new pfloat[m];
  h[0] = -0.5*n; // -n/2
  for (int i=1;i<m;i++) { h[i] = 0.0; }
  pfloat* b = NULL;

  pwork* mywork = ECOS_setup(n, m, p, l, ncones, q, e, Gpr, Gjc, Gir, Apr, Ajc, Air, C, h, b);
  if (mywork != NULL) {
    int exitCode = ECOS_solve(mywork);

    for (int i=0;i<n;i++) { cout<<mywork->best_x[i]<<" "; } //value of the variale
    cout<<exitCode<<"\n";
    ECOS_cleanup(mywork, 0);
  } else {
    cout<<"Could not initialize ECOS_setup\n";
  }
  return 0;
}
