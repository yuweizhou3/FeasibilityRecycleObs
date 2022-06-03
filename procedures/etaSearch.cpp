#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cmath>

#define alpha  0.05
#define k        1
#define n        20
#define C        1
#define high    2.0
#define low     0.0


double check(double,double);

FILE *outfile;

int main()
{
    outfile = NULL;
    outfile = fopen("eta.out","a");


   double h,l, half;
   double value_h, value_l;

   fprintf(outfile, "\n\n n: %d \t k : %d\n",n,k);

   for(int i= 1; i <= C; i++){

   h = high;
   l = low;

   value_h = check(h,i);
   value_l = check(l,i);

   if (value_h < 0 && value_l > 0){
       while( fabs( h - l) > 0.0001 ){
           half = (h+l)/2.0;
           if( check(half,i) > 0 )   l = half;
           else                             h = half;
       }

       fprintf(outfile, " c: %d \t",i);
      fprintf(outfile, "eta is between %.4f and %.4f\n", l,h);

   }

   else  fprintf(outfile, "Error! Error! Error!\n");

   }
   fclose(outfile);

}


double check(double eta, double c)
{
    double temp = 0, factor = 1;

   for(int l = 1; l <= c; l++){
      if( l == c)  factor = 0.5;
      temp += pow(-1, l+1)*factor*pow(1.0 + 2*eta*(2*c-l)*l/c, (double)-(n-1)/2);
   }
   //printf("%.10f\n", pow(1-alpha, 1/k));

   //return(temp - alpha/(k-1));
   //return (temp - alpha/k);
   return (temp - (0.05));
}


