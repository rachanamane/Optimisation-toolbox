#include<iostream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

#include "ecos.h"

pfloat getNextDouble(ifstream &infile, bool isEndOfLine) {
  char cNum[100];
  if (isEndOfLine) {
    infile.getline(cNum, 256);
  } else {
    infile.getline(cNum, 256, ',');
  }
  return atof(cNum);
}

void readArray(ifstream &infile, pfloat *array, int size) {
  for (int i=0;i<size;i++) {
    array[i] = getNextDouble(infile, (i==(size-1)));
  }
}

int main () {
  ifstream infile;
  infile.open("sample.txt", ifstream::in);
  if (!infile.is_open()) {
    cout<<"Could not open file"<<endl;
    return 0;
  }
  int n = getNextDouble(infile, false);
  int m = getNextDouble(infile, true);
  // minimize: c.x
  // subject to: s <= Ax <= b
  // l <= x <= u
  pfloat *C = new pfloat[n]; readArray(infile, C, n);
  pfloat *s = new pfloat[m]; readArray(infile, s, m);
  pfloat *b = new pfloat[m]; readArray(infile, b, m);
  pfloat *l = new pfloat[n]; readArray(infile, l, n);
  pfloat *u = new pfloat[n]; readArray(infile, u, n);
  pfloat **A = new pfloat*[m];
  for (int i=0;i<m;i++) {
    A[i] = new pfloat[n];
    readArray(infile, A[i], n);
  }
  infile.close();
  cout<<"\nFile read successfully\n";

  // Inequalities: 
  // Ax <= b: m rows
  // -Ax <= -s: m rows
  // x <= u: n rows
  // -x <= -l: n rows
  int inequalityCount = 2*(m+n);

  // G is size (2*m+2n)xn
  // https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_.28CSC_or_CCS.29
  // Creating a temporary vector since we don't know the size of Gpr and Gir to create an array.
  vector<pfloat> GprVector;
  idxint* Gjc = new idxint[n+1];
  vector<pfloat> GirVector;
  pfloat* h = new pfloat[inequalityCount];
  int nonZeroCount = 0;
  Gjc[0] = 0;
  for (int j=0;j<n;j++) {  // Iterating over columns to move from 0th to nth column
    for (int i=0;i<inequalityCount;i++) {
      pfloat gVal;
      if (i<m) {
	    // Ax<=b
        gVal = A[i][j];
        if (j==0) { h[i] = b[i]; }
      } else if (i<2*m) {
        // -Ax <= -s
        gVal = -A[i-m][j];
        if (j==0) { h[i] = -s[i-m]; }
      } else if (i<(2*m+n)) {
        // x <= u (LHS is diagonal matrix of size n)
        gVal = (i-2*m == j) ? 1.0 : 0.0;
        if (j==0) { h[i] = u[i-2*m]; }
      } else {
        // -x <= -l  (LHS is diagonal matrix of size n)
        gVal = (i-2*m-n == j) ? -1.0 : 0.0;
        if (j==0) { h[i] = -l[i-2*m-n]; }
      }
      if (gVal != 0) {  // Float comparision??????
        nonZeroCount++;
        GprVector.push_back(gVal);
        GirVector.push_back(i);
      }
    }
    Gjc[j+1] = nonZeroCount;
  }
  pfloat* Gpr = new pfloat[GprVector.size()];
  idxint* Gir = new idxint[GirVector.size()];
  for (int i=0;i<(int) GprVector.size();i++) {
    Gpr[i] = GprVector.at(i);
    Gir[i] = GirVector.at(i);
  }
  pwork* mywork = ECOS_setup(n, inequalityCount, 0, inequalityCount, 0, new idxint[0], 0, Gpr, Gjc, Gir, NULL, NULL, NULL, C, h, NULL);
  if (mywork != NULL) {
    cout<<"ECOS setup successful. Solving now...";
    int exitCode = ECOS_solve(mywork);
    cout<<"Status returned by solver: "<<exitCode<<endl;
    for (int i=0;i<n;i++) { cout<<"x"<<i+1<<": "<<mywork->best_x[i]<<endl; }
    cout<<endl;
    ECOS_cleanup(mywork, 0);
  } else {
    cout<<"Could not initialize ECOS_setup\n";
  }
  return 0;
}
