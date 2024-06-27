
//
//  Created by Osakana on 2022/05/16.
//  v2: 2022/5/22.
//

#include <stdio.h>
#include <math.h>
#include "Legendre-zero-points/GaussLegendreIntegral.h"

#define NORMAL 10. // normalization factor (or, max value of integral variable)

double integrand(double x)
{
	return NORMAL*exp(-NORMAL*NORMAL*x*x);
}

double integrand2(double x[])
{
	return pow(NORMAL,2)*exp(-pow(NORMAL,2)*(x[0]*x[0]+x[1]*x[1]));
}

double integrand4(double x[])
{
	return pow(NORMAL,4)*exp(
				   -pow(NORMAL,2)*(
						 x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]
						 )
				   );
}







int main(int argc, const char * argv[]) {
	
	
	double x[4];
	double I=0.;
	
	printf("%e\n",x);
	I = GLintegral(integrand,8);// would return sqrt(pi)
	printf("result = %e,\t",I);
	printf("error = %e\n",I*I-M_PI);
	
	I = GLintegral4(integrand4,x,0,1,2,3,6); // would return pi^2
	printf("result = %e,\t",I);
	printf("error = %e\n",I-M_PI*M_PI);
	
	return 0;
}
