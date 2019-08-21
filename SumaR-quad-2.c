// Suma de Riemann de sigma BII

#include <stdio.h>
#include <math.h>
#include <quadmath.h>
#include "pralq.h"


#define MASA 1.
#define CARGA 1.

#define OMEGA 100.


/* #define PARTP2 */

#define PARTP2 1000.
#define PARTTHETA1 200.
#define PARTP1 100.

	FILE *fptr2;


__float128 bucle1(__float128 p1, __float128 theta1, __float128 nveces)
{
	__float128 number = 0.;

if(det(MASA,OMEGA,p1,theta1)>0.)
{

	/* Partición de p2*/
	__float128 deltap2=(p2b(MASA,OMEGA,p1,theta1)-p2a(MASA,OMEGA,p1,theta1))/PARTP2;
        __float128 p2 = p2a(MASA,OMEGA,p1,theta1)*1.0000001+ nveces*((deltap2*PARTP2)/100.); /* En p2a hay problema, hace falta correr un poquito el inicio de la suma*/


        do
        {
            if(deltap2>0.)
	    {
	        number += sigma(CARGA,MASA,OMEGA,p1,theta1,p2)*deltap2;
            }
	    p2 += deltap2;
        }
        while(p2 <= p2b(MASA,OMEGA,p1,theta1)-((99.-nveces)*(deltap2*PARTP2)/100.));
}
	return number;
}

__float128 bucle2(__float128 p1,__float128 nveces)
{
	/* Particion de theta1 */

	__float128 number1 = 0., theta1=0., bb = 0., deltatheta1=theta1max(MASA,OMEGA)/PARTTHETA1;


	do
	{
	    bb=bucle1(p1,theta1,nveces);
	    number1 += bb*deltatheta1;
	    printf("\n theta1=%Qf \t bb=%0.20Qf \n",theta1,bb);
       	    fprintf(fptr2,"%Qf \t %0.20Qf \n",theta1,bb);
	    theta1 += deltatheta1;
	 }
	 while(theta1 <= theta1max(MASA,OMEGA));

	return number1;
}

__float128 bucle3(__float128 nveces)
{

	/* Particion de p1 */ 

	__float128 number3 = 0., cc = 0., deltap1= (p1max(MASA,OMEGA)-p1min(MASA,OMEGA))/PARTP1 , p1 = p1min(MASA,OMEGA) ;
	
	//FILE *fptr;
	//fptr = fopen("datos-sr.dat","w");
	do
	{
		cc = bucle2(p1,nveces);

                if(isnan(cc)) {cc=0; }

		number3 += cc*deltap1;
		//printf("\n p1=%Qf \t cc=%0.20Qf \n",p1,cc);
		//fprintf(fptr,"%Qf \t %0.20Qf \n",p1,cc);
		p1 += deltap1;
	}
	while(p1 <= p1max(MASA,OMEGA)*1.0001);
	//fclose(fptr);
	return number3;   

}

int main()
{

    //printf("\n Suma de Riemann = %Qf \n", bucle3(1));


__float128 number4 = 0. ,dd=0. , nveces=0.;

	FILE *fptr;
	fptr = fopen("datos-partes.dat","w");
	fptr2 = fopen("datos-theta1.dat","w");


	do
	{
		dd = bucle3(nveces);

		number4 += dd;
		printf("\n n=%Qf \t cc=%0.20Qf \n",nveces,number4);
		fprintf(fptr,"%Qf \t %0.20Qf \n",nveces,number4);
		nveces += 1.;
		}
	while(nveces < 1.);
	fclose(fptr);
	fclose(fptr2);
 
    return 0;
}
