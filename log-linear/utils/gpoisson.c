#include <R.h>
#include <Rmath.h>

double fullloglik(int *X, int *Y, int size_x, int size_y, double *parameters)
{
    int i;
    double result = 1.0;
    
    for (i = 0; i < size_x; i++)
        result *= dpois(X[i], parameters[i], FALSE);
    
    for (i = 0; i < size_y; i++)
        result *= dpois(Y[i], parameters[size_x + i], FALSE);
    return(result);
}

double lik(int *X, int *Y, int size_x, int size_y, int r, double *parameters)
{
    int i, j, s,e, nr, ymax, NX[size_x], NY[size_y];
    double result,val;
    
    if (r == size_y) {
        result = fullloglik(X, Y, size_x, size_y, parameters);
    } else {
        for (s=0,j=0,i=j+1;s<r;i++) {
          s++;
          if (i==size_x) {
            j++;i=j+1;
          }
        } 
           
        nr = r + 1;
        ymax = fmin(X[i], X[j]);
        result = 0;
        for (e = 0; e <= ymax; e++) {
            
            memcpy(NX, X, size_x * sizeof(int));
            NX[i] =  NX[i] - e;
            NX[j] =  NX[j] - e;
            
            memcpy(NY, Y, size_y * sizeof(int));
            NY[nr-1] = e;
            val=lik(NX, NY, size_x, size_y, nr, parameters);
            result += val;
       }
    }
    return(result);
}

void gpoisson(int *X, int *Y, int *size_x, int *size_y, double *parameters, double *ans)
{
    *ans = log(lik(X, Y, *size_x, *size_y, 0, parameters));
}
