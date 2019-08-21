// Suma de Riemann de sigma

#include <stdio.h>
#include <math.h>
#include <quadmath.h>
#include "pralq.h"

__float128 bucle1(__float128 p1, __float128 theta1)
{
	__float128 number = 0.;

if(det(1.,100.,p1,theta1)>0.)
{

	/* Partición de p2*/
	__float128 deltap2=(p2b(1.,100.,p1,theta1)-p2a(1.,100.,p1,theta1))/100.;
    __float128 p2 = p2a(1.,100.,p1,theta1)+deltap2/100000.; /* En p2a hay problema, hace falta correr un poquito el inicio de la suma*/


        do
        {
            if(deltap2>0.)
	    {
	        number += sigma(1.,1.,100.,p1,theta1,p2)*deltap2;
		//printf("p2=%f \t" ,p2,theta1,p2);
		//number += deltap2;
		//printf("suma=%f \n" ,number);
            }
	    p2 += deltap2;
        }
        while(p2 <= p2b(1.,100.,p1,theta1));
}
	return number;
}

__float128 bucle2(__float128 p1)
{
	/* Particion de theta1 */

	__float128 number1 = 0., theta1=0., bb = 0., deltatheta1=theta1max(1.,100.)/100.;
	do
	{
	    bb=bucle1(p1,theta1);
	    number1 += bb*deltatheta1;
	    theta1 += deltatheta1;
	 }
	 while(theta1 <= theta1max(1.,100.));
	return number1;
}

__float128 bucle3(void)
{

	/* Particion de p1 */ 

	__float128 number3 = 0., p1 = p1min(1.,100.) , cc = 0., deltap1= (p1max(1.,100.)-p1min(1.,100.))/50.;
	
	FILE *fptr;
	fptr = fopen("datos-sr.dat","w");
	do
	{
		cc = bucle2(p1);

		number3 += cc*deltap1;
		printf("\n p1=%Qf \t cc=%0.20Qf \n",p1,cc);
		fprintf(fptr,"%Qf \t %0.20Qf \n",p1,cc);
		p1 += deltap1;
	}
	while(p1 <= p1max(1.,100.)+p1max(1.,100.)/100000.);
	fclose(fptr);
	return number3;   

}

int main()
{

//printf("bucle1 = %f \n", bucle1(1,0.1));

//printf("\n \n Suma de Riemann = %lf \n", bucle2(4));

    printf("\n Suma de Riemann = %Qf \n", bucle3());
 
    return 0;
}
