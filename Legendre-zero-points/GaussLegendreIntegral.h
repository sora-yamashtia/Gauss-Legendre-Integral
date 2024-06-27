//
//  GaussLegendreIntegral.h
//  #include "Legendre-zero-points/GaussLegendreIntegral.h"
//
//  Created by osakana on 2023/09/15.
//




#ifndef GaussLegendreIntegral_h
#define GaussLegendreIntegral_h


double LegendreWeight(int n, double x);
double GLintegral(double (*f)(double),int n);
double GLintegral2(double (*f)(double x[]), double x[],
				   int a, int b, int n);
double GLintegral3(double (*f)(double x[]), double x[],
				   int a, int b, int c, int n);
double GLintegral4(double (*f)(double x[]), double x[],
				   int a, int b, int c, int d, int n);

double LegendreWeight(int n, double x)
{
	double p1=1.;// P_n
	double p2=0.;// P_{n-1}
	double p3=0.;// P_{n-2}
	
	// using P_n = 1/n[ (2n-1)xP_{n-1} -(n-1)P_{n-2} ]
	for(int i=1;i<=n;i++)
	{
		p3=p2;
		p2=p1;
		p1 = ((2.*i-1)*x*p2-(i-1.)*p3)/i;
	}
	
	// using P_n'(x) = n/(x^2-1) * [xP_n - P_{n-1} ]
	p1 = n*(x*p1-p2)/(x*x-1.);
	
	// get weight
	p1 = 2./(1.-x*x)/(p1*p1);
	
	return p1;
}

double GLintegral(double (*f)(double),int n)
{
	int precision = round(pow(2,n));
	printf("precision = %d\n",precision);
	
	int i=0;
	double variable;
	double w;//weight
	char filename[128];
	sprintf(filename,"Legendre-zero-points/Legendrezero-%d.dat",n);
	printf("LOADING %s ...\n\n\n",filename);
	
	FILE *file;
	
	file=fopen(filename,"r");
	
	
	double integral=0.;
	i=0;
	while(fscanf(file, "%lf",&variable) != EOF)
	{
		printf("x = %e\n",variable);
		
		w = LegendreWeight(precision,variable);
		integral = integral + f(variable)*w;
		
	}
	
	printf("integral_end %e\n",M_PI - integral*integral);
	return integral;
}


double GLintegral2(double (*f)(double x[]), double x[], int a, int b, int n)
{
	int precision = round(pow(2,n));
	printf("precision = %d\n",precision);
	
	int i=0;
	int j=0;
	double variable[precision];
	double w[precision];
	char filename[128];
	sprintf(filename,"Legendre-zero-points/Legendrezero-%d.dat",n);
	printf("LOADING %s ...\n\n\n",filename);
	
	FILE *file;
	
	file=fopen(filename,"r");
	
	
	double integral=0.;
	i=0;
	while(fscanf(file, "%lf",&variable[i]) != EOF)
	{
		w[i] = LegendreWeight(precision,variable[i]);
		printf("x = %+e, weight=%+e\n",variable[i],w[i]);
		
		i++;
	}
	
	integral=0;
	for (i=0; i<precision; i++)
	{
		for (j=0; j<precision; j++)
		{
			x[a] = variable[i];
			x[b] = variable[j];
			integral = integral + w[i]*w[j]*f(x);
			
			//printf("integral = %e\n",integral);
		}
	}
	
	
	
	printf("integral_end %e\n",M_PI - integral);
	return integral;
}

double GLintegral3(double (*f)(double x[]), double x[],
				   int a, int b, int c, int n)
{
	int precision = round(pow(2,n));
	printf("precision = %d\n",precision);
	
	int i=0;
	int j=0;
	int k=0;
	double variable[precision];
	double w[precision];
	char filename[128];
	sprintf(filename,"Legendre-zero-points/Legendrezero-%d.dat",n);
	printf("LOADING %s ...\n\n\n",filename);
	
	FILE *file;
	
	file=fopen(filename,"r");
	
	
	double integral=0.;
	i=0;
	while(fscanf(file, "%lf",&variable[i]) != EOF)
	{
		w[i] = LegendreWeight(precision,variable[i]);
		printf("x = %+e, weight=%+e\n",variable[i],w[i]);
		
		i++;
	}
	
	integral=0;
	for (i=0; i<precision; i++)
	{
		for (j=0; j<precision; j++)
		{
			for (k=0;k<precision; k++)
			{
				x[a] = variable[i];
				x[b] = variable[j];
				x[c] = variable[k];
				integral = integral + w[i]*w[j]*w[k]*f(x);
				
				//printf("integral = %e\n",integral);
				
			}
		}
	}
	
	
	
	printf("integral_end %e\n",M_PI - integral);
	return integral;
}




double GLintegral4(double (*f)(double x[]), double x[],
				   int a, int b, int c, int d, int n)
{
	int precision = round(pow(2,n));
	printf("precision = %d\n",precision);
	
	int i=0;
	int j=0;
	int k=0;
	int l=0;
	double variable[precision];
	double w[precision];
	char filename[128];
	sprintf(filename,"Legendre-zero-points/Legendrezero-%d.dat",n);
	printf("LOADING %s ...\n\n\n",filename);
	
	FILE *file;
	
	file=fopen(filename,"r");
	
	
	double integral=0.;
	i=0;
	while(fscanf(file, "%lf",&variable[i]) != EOF)
	{
		w[i] = LegendreWeight(precision,variable[i]);
		printf("x = %+e, weight=%+e\n",variable[i],w[i]);
		
		i++;
	}
	
	integral=0;
	for (i=0; i<precision; i++)
	{
		for (j=0; j<precision; j++)
		{
			for (k=0;k<precision; k++)
			{
				for (l=0; l<precision; l++)
				{
					x[a] = variable[i];
					x[b] = variable[j];
					x[c] = variable[k];
					x[d] = variable[l];
					integral = integral + w[i]*w[j]*w[k]*w[l]*f(x);
					
					//printf("integral = %e\n",integral);
					
				}
				
			}
		}
	}
	
	
	
	printf("integral_end %e\n",M_PI*M_PI - integral);
	return integral;
	
}

#endif /* GaussLegendreIntegral_h */
