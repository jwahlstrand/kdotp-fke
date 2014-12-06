#ifndef TENSORS_H
#define TENSORS_H

#define X 0
#define Y 1
#define Z 2

typedef struct {
    double diag[3];
    complex double offdiag[3];
} Tensor2;

typedef struct {
    complex double data[27];
} Tensor3;

#define incr3(t,i,j,k,val) (t.data[i*9+j*3+k]+=val)

typedef struct {
    complex double data[81];
} Tensor4;

#define incr4(t,i,j,k,l,val) (t.data[i*27+j*9+k*3+l]+=val)
#define get4(t,i,j,k,l) (t.data[i*27+j*9+k*3+l])

#endif

