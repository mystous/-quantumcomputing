#include <iostream>
#include <cmath>

using namespace  std;

void matrixvectormul(double *matrix, int matrixM, int matrixN, double *vector, int vectorDim, double *result, int resultDim)
{
    int i, j;

    if ( resultDim < vectorDim )
        return;

    if ( matrixN != vectorDim )
        return;

    for ( i = 0 ; i < resultDim ; ++i )
        result[i] = 0.0;

    for ( i = 0 ; i < matrixM ; ++i )
    {
        for ( j = 0 ; j < matrixN ; ++j )
        {
            result[i] += (matrix[ i*matrixN + j ]  * vector[j]);
        }
    }
}

void vectorSwap(double *vectorSource, double *vectorTarget, int vectorDim)
{
    int i;

    for ( i = 0 ; i < vectorDim ; ++i )
        vectorTarget[i] = vectorSource[i];
}

void tensorProduct(double *tensor1, int tensor1DimX, int tensor1DimY, double *tensor2, int tensor2DimX, int tensor2DimY,  double *result, int resultDimX, int resultDimY)
{
    int t1x, t1y, t2x, t2y;

    for ( t1x = 0 ; t1x < tensor1DimX ; ++t1x )
    {
        for ( t1y = 0 ; t1y < tensor1DimY ; ++t1y )
        {
    
            for ( t2x = 0 ; t2x < tensor2DimX ; ++t2x )
            {
                for ( t2y = 0 ; t2y < tensor2DimY ; ++t2y )
                {
                    result[ ( ((t1y*tensor2DimY) + t2y) * (tensor1DimX*tensor2DimX) ) + (t1x*tensor2DimX + t2x) ] = tensor1[t1y*tensor1DimX+t1x] * tensor2[t2y*tensor2DimX+t2x];
                }
            }
        }
    }

}

int main() {
    int i;
    double q0[2], q1[2], qdumy[4], q01[4];
    double entq2[4], entq2p[4];
    double *hadamard = new double[2*2];
    double *hadamardp = new double[4*4];
    double *controlNotGate = new double[4*4];

    q0[0] = q1[0] = 1.0;
    q0[1] = q1[1] = 0.0;
    qdumy[0] = qdumy[1] = 0.0;

    hadamard[0] = hadamard[1] = hadamard[2] = 1/sqrt(2);
    hadamard[3]= -1/sqrt(2);

    for ( i = 0 ; i < (4 * 4); ++i )
        controlNotGate[i] = 0.0;
    controlNotGate[0] = controlNotGate[5] = 1;
    controlNotGate[11] = controlNotGate[14] = 1;

    tensorProduct(hadamard, 2, 2, hadamard, 2, 2, hadamardp, 4, 4);
    tensorProduct(q0, 1, 2, q1, 1, 2, q01, 1, 4);

    matrixvectormul(hadamardp, 4, 4, q01, 4, qdumy, 4);
    vectorSwap(qdumy, q01, 4);
 
    matrixvectormul(controlNotGate, 4, 4, q01, 4, entq2, 4);

    for ( i = 0 ; i < 4 ; ++i )
        entq2p[i] = pow(entq2[i], 2);
    printf("Pr(|00>) = %f, Pr(|01>) = %f, Pr(|10>) = %f, Pr(|11>) = %f\n", entq2p[0], entq2p[1], entq2p[2], entq2p[3]);
 
    delete hadamard;
    delete hadamardp;
    delete controlNotGate;
} 
