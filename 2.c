#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double num_comput_integral_l_re(double left_boundary_a,
double right_boundary_b, unsigned int intervals);
double num_comput_integral_r_re(double left_boundary_a,
double right_boundary_b, unsigned int intervals);
double num_comput_integral_Simps (double left_boundary_a,
double right_boundary_b, unsigned int intervals);
double num_comput_integral_are_of_the_trapecium (double left_boundary_a,
double right_boundary_b, unsigned int intervals);
double integrand_expression( double x );
double max_error = 0.00001;
    int main()
{
 double left_boundary_a=0, right_boundary_b=1;
 double measurement_error=0, I1=0, I2=0;
 int intervals, var, i,errors;
 double integral_s=0;



   while(1)
 {
 printf("\n\tEnter the left boundary X1 =");
 scanf("%lf", &left_boundary_a);
 printf("\n\tEnter the right boundary X2 =");
 scanf("%lf", &right_boundary_b);
 do
 {
 printf("\tEnter the number of partition intervals (N>0)\nN=");
 scanf("%u", &intervals);
 }
    while(intervals <= 0);
    do
{
 printf("\nChoose the method of calculating:\n");
 printf("\t1. By Left Rectangles :\n");
 printf("\t2. By Right Rectangles:\n");
 printf("\t3. By ntegral_Simps's method (parabola method):\n");
 printf("\t4. By are of the trapecium:\n ");
 scanf("%u", &var);
 if (var!=1 && var!=2 &&var!=3 && var!=4)
 printf("\nYou are mistaken\n");
 }
 while (var!=1 && var!=2 && var!=3 && var!=4);
 system("cls");
 switch(var)
 {
 case 1:
 {
   integral_s= num_comput_integral_l_re( left_boundary_a, right_boundary_b,intervals);
   printf("\n\n\t======Right Rectangles method======\n");
   printf("\n\ta = %.2lf \n\tb = %.2lf \n\tIntegral = %.8lf \n\tN = %d",
   left_boundary_a, right_boundary_b,integral_s,intervals);
   while(1){
      intervals = intervals+2;
      float I1 = num_comput_integral_l_re(left_boundary_a,right_boundary_b,intervals);
      float I2 = num_comput_integral_l_re(left_boundary_a,right_boundary_b,intervals+2);
      if(fabs(I2-I1)<max_error){break;}
   }
   float n_integral_s = num_comput_integral_are_of_the_trapecium(left_boundary_a,right_boundary_b,intervals);
   float I1 = integral_s - n_integral_s;
   printf("\n\tN1:%d\n\tIntegral:%8lf\n\tI1:%8lf", intervals,n_integral_s,I1);
   }
 break;
 case 2:
 {
   integral_s= num_comput_integral_r_re( left_boundary_a, right_boundary_b,intervals);
   printf("\n\n\t======Right Rectangles method======\n");
   printf("\n\ta = %.2lf \n\tb = %.2lf \n\tIntegral = %.8lf \n\tN = %d",
   left_boundary_a, right_boundary_b,integral_s,intervals);
   while(1){
      intervals = intervals+2;
      float I1 = num_comput_integral_r_re(left_boundary_a,right_boundary_b,intervals);
      float I2 = num_comput_integral_r_re(left_boundary_a,right_boundary_b,intervals+2);
      if(fabs(I2-I1)<max_error){break;}
   }
   float n_integral_s = num_comput_integral_r_re(left_boundary_a,right_boundary_b,intervals);
   float I1 = integral_s - n_integral_s;
   printf("\n\tN1:%d\n\tIntegral:%8lf\n\tI1:%8lf", intervals,n_integral_s,I1);
 }
 break;
 case 3:
 {
   integral_s= num_comput_integral_Simps( left_boundary_a, right_boundary_b,intervals);
   printf("\n\n\t======Right Rectangles method======\n");
   printf("\n\ta = %.2lf \n\tb = %.2lf \n\tIntegral = %.8lf \n\tN = %d",
   left_boundary_a, right_boundary_b,integral_s,intervals);
   while(1){
      intervals = intervals+2;
      float I1 = num_comput_integral_Simps(left_boundary_a,right_boundary_b,intervals);
      float I2 = num_comput_integral_Simps(left_boundary_a,right_boundary_b,intervals+2);
      if(fabs(I2-I1)<max_error){break;}
   }
   float n_integral_s = num_comput_integral_Simps(left_boundary_a,right_boundary_b,intervals);
   float I1 = integral_s - n_integral_s;
   printf("\n\tN1:%d\n\tIntegral:%8lf\n\tI1:%8lf", intervals,n_integral_s,I1);
 }
 break;
 case 4:
 { 
  integral_s= num_comput_integral_are_of_the_trapecium( left_boundary_a, right_boundary_b,intervals);
   printf("\n\n\t======Right Rectangles method======\n");
   printf("\n\ta = %.2lf \n\tb = %.2lf \n\tIntegral = %.8lf \n\tN = %d",
   left_boundary_a, right_boundary_b,integral_s,intervals);
   while(1){
      intervals = intervals+2;
      float I1 = num_comput_integral_are_of_the_trapecium(left_boundary_a,right_boundary_b,intervals);
      float I2 = num_comput_integral_are_of_the_trapecium(left_boundary_a,right_boundary_b,intervals+2);
      if(fabs(I2-I1)<max_error){break;}
   }
   float n_integral_s = num_comput_integral_are_of_the_trapecium(left_boundary_a,right_boundary_b,intervals);
   float I1 = integral_s - n_integral_s;
   printf("\n\tN1:%d\n\tIntegral:%8lf\n\tI1:%8lf", intervals,n_integral_s,I1);
 }
break;
   }
   break;
}
 return 0;
}

double num_comput_integral_l_re(double left_boundary_a,double right_boundary_b, unsigned int intervals)
{
 double integral_s = 0, x = 0, h;
  unsigned int i;
  h = (right_boundary_b - left_boundary_a) / intervals;
  x = left_boundary_a;
  for (i = 0; i <= (intervals - 1); i++)
  {
   integral_s += integrand_expression(x);
   x += h;
  }
  return integral_s * h;
}
double num_comput_integral_r_re(double left_boundary_a,double right_boundary_b, unsigned int intervals)
{
   double integral_s = 0, x = 0, h;
   unsigned int i;
   h = (right_boundary_b - left_boundary_a) / intervals;
   x = left_boundary_a + h;
   for (i = 0; i <= intervals; i++)
   {
    integral_s += integrand_expression(x);
    x += h;
   }
   return integral_s * h;
   
}
double num_comput_integral_Simps (double left_boundary_a,double right_boundary_b, unsigned int intervals)
{
    double integral_s=0, h=0, x, sum1=0,sum2=0;
    unsigned int i;
    h = (right_boundary_b - left_boundary_a) / intervals;
    x=left_boundary_a + i*h;
    for(i=1;i<(intervals-1);i++)
    {
    if(i%2==0)
      sum1 += integrand_expression(left_boundary_a + h * i);
    else
      sum2 += integrand_expression(left_boundary_a+ h * i);
  }
  return (h/3)*(integrand_expression(left_boundary_a) + integrand_expression(right_boundary_b) + 4 * sum1 + 2 * sum2);
    
}
double num_comput_integral_are_of_the_trapecium(double left_boundary_a,double right_boundary_b,unsigned intervals)
{
  double integral_s = 0, x = 0, h;
  unsigned int i;
  h = (right_boundary_b - left_boundary_a) / intervals;
  x = left_boundary_a + h;
  for (i = 0; i <= (intervals - 1); i++)
  {
    integral_s += (integrand_expression(x) + integrand_expression(x + h)) / 2;
    x += h;
  }
  return integral_s * h;
}
double integrand_expression( double x )
{
   return pow(x,2)*sqrt(1+pow(x,3));
}
