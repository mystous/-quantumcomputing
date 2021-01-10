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

void tensorProduct(double *vector1, int vector1Dim, double *vector2, int vector2Dim, double *result, int resultDim)
{
    int i, j;

    for ( i = 0 ; i < vector1Dim ; ++i )
    {
        for( j = 0 ; j < vector2Dim ; ++j )
        {
            if( (i*vector2Dim + j) >= resultDim )
                break;
            result[ i*2 + j ] = vector1[i] * vector2[j];
        }
    }
}

int main() {
    int i;
    double q0[2], q1[2], qdumy[2];
    double entq2[4], entq2p[4];
    double *hadamard = new double[2*2];
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

    matrixvectormul(hadamard, 2, 2, q0, 2, qdumy, 2);
    vectorSwap(qdumy, q0, 2);

    matrixvectormul(hadamard, 2, 2, q1, 2, qdumy, 2);
    vectorSwap(qdumy, q1, 2);

    tensorProduct(q0, 2, q1, 2, entq2, 4);

    matrixvectormul(controlNotGate, 4, 4, entq2, 4, entq2p, 4);

    for ( i = 0 ; i < 4 ; ++i )
        entq2p[i] = pow(entq2p[i], 2);
    printf("Pr(|00>) = %f, Pr(|01>) = %f, Pr(|10>) = %f, Pr(|11>) = %f\n", entq2p[0], entq2p[1], entq2p[2], entq2p[3]);
 
    delete hadamard;
    delete controlNotGate;
} 
