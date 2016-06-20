#include "computeRoeSolver.h"
#include "macroDefinition.h"

#include <math.h>
#include <stdio.h>

void computeRoeSolver(double hl, double hr, double ul, double ur, double vl, double vr,
					 int g_n, int f_n, double *F, double *amax)
{
	double duml, dumr, cl, cr;
	double uperp, h_hat, u_hat, v_hat, c_hat;
	double dh, du, dv, dupar, duperp;
	double res[3]={0.0};
	// double multRA[3][3] = {{0.0}};
	double da1, da3, a1, a2, a3, al1, al3, ar1, ar3;
	double uperpl, uperpr;
	duml = sqrt(hl);
	//printf("duml: %lf\n", duml);
	dumr = sqrt(hr);
	//printf("dumr: %lf\n", dumr);
	cl = sqrt(gravity * hl);
	cr = sqrt(gravity * hr);
	h_hat = duml*dumr;
	u_hat = (duml*ul + dumr*ur) / (dumr + duml);
	v_hat = (duml*vl + dumr*vr) / (duml + dumr);
	c_hat = sqrt(0.5*gravity*(hl + hr));
	uperp = u_hat*f_n + v_hat*g_n;
	dh = hr - hl;
	du = ur - ul;
	dv = vr - vl;

	dupar = -du * g_n + dv * f_n;
	duperp = du * f_n + dv * g_n;
			// printf("dupar: %lf\n", dupar);

	/* R matrix [1 0 1;
	uhat-chat*cn -sn uhat+chat*cn;
    vhat-chat*sn cn vhat+chat*sn]*/
	double R[3][3] = {{1.0, 0.0, 1.0}, {u_hat - c_hat*f_n, -g_n, u_hat + c_hat*f_n }, {v_hat - c_hat*g_n, f_n, v_hat + c_hat*g_n}};
    // double R[3][3] = {{1.0, 0.0, 0.0}, {0., 1., 0.}, {0., 0., 1.}};
    /* R^-1 * dQ = dW*/
	double dW[3] = {0.5 * (dh - h_hat*duperp/c_hat),
		  h_hat * dupar,
		  0.5 * (dh + h_hat*duperp/c_hat)};
		  // printf("dW[0]: %lf\n", dW[0]);
    /*Critical flow fix*/
	uperpl = ul*f_n + vl*g_n;
	uperpr = ur*f_n + vr*g_n;
    al1 = uperpl-cl;
	al3 = uperpl+cl;
	ar1 = uperpr-cr;
	ar3 = uperpr+cr;
    da1 = fmax(0.0, 2*(ar1 - al1));
    da3 = fmax(0.0, 2*(ar3 - al3));
    a1 = fabs(uperp - c_hat);
    a2 = fabs(uperp);
    a3 = fabs(uperp + c_hat);

    if (a1<da1)
    {
    	a1 = 0.5*(a1*a1/da1 + da1);
    }
    if (a3<da3)
    {
    	a3 = 0.5*(a3*a3/da3 + da3);
    }

    /*Compute flux*/
    double A[3][3] = {{a1, 0.0, 0.0},
         {0.0, a2, 0.0},
   	     {0.0, 0.0, a3}};
   	//FL[3] = { uperpl*hl, ul*uperpl*hl + 0.5*gravity*hl*hl*f_n, vl*uperpl*hl + 0.5*gravity*hl*hl*g_n};
   	double FL[3] = { uperpl*hl, ul*uperpl*hl + 0.5*gravity*hl*hl*f_n, vl*uperpl*hl + 0.5*gravity*hl*hl*g_n};
   	double FR[3] = { uperpr*hr, ur*uperpr*hr + 0.5*gravity*hr*hr*f_n, vr*uperpr*hr + 0.5*gravity*hr*hr*g_n};
 	// printf("f_n:%d\n", f_n);
 	// printf("hr:%lf\n", hr);
 	// printf("uperpr:%lf\n", uperpr);
 	// printf("ur*uperpr*hl:%lf\n", ur*uperpr*hl);
 	// for (int i = 0; i < 3; ++i)
 	// {
 	// 	printf("FR[%d]:%lf\n", i, FR[i]);
 	// }
 	// printf("\n");
 //   		for (int i = 0; i < 3; ++i)
	// {	
	// 	for (int j = 0; j < 3; ++j)
	// 	{
	// 		printf("A[%d][%d]: %lf\t",i, j, A[i][j]);
	// 	}
	// 	printf("\n");
	// }
	// 	for (int i = 0; i < 3; ++i)
	// {
	// 	for (int k = 0; k < 3; ++k)
	// 	{
	// 		for (int j = 0; j < 3; ++j)
	// 		{
	// 			multRA[i][j] += R[i][k] * A[k][j];
	// 		}
	// 	}
	// }
	// for (int i = 0; i < 3; ++i)
	// {	
	// 	for (int j = 0; j < 3; ++j)
	// 	{
	// 		printf("multRA[%d][%d]: %lf\t",i, j, multRA[i][j]);
	// 	}
	// 	printf("\n");
	// }
	// for (int i = 0; i < 3; ++i)
	// {
	// 	for (int j = 0; j < 3; ++j)
	// 	{	
	// 		res[i] += multRA[i][j] * dW[j];
	// 	}
	// }
   	computeRes(R, A, dW, res);
   	// printf("res[0]: %lf\n", res[0]);
   	for (int i = 0; i < 3; ++i)
   	{
   	   	F[i] = 0.5 * (FL[i] + FR[i] - res[i]);
   	}
   	//printf("c_hat: %lf\n", c_hat);
   	*amax = c_hat + fabs(uperp);

}

void computeRes(double R[3][3], double A[3][3], double *dW, double *res){
	double multRA[3][3] = {{0.0}};
	for (int i = 0; i < 3; ++i)
	{
		for (int k = 0; k < 3; ++k)
		{
			for (int j = 0; j < 3; ++j)
			{
				multRA[i][j] += R[i][k] * A[k][j];
			}
		}
	}
	// for (int i = 0; i < 3; ++i)
	// {	
	// 	for (int j = 0; j < 3; ++j)
	// 	{
	// 		printf("multRA[%d][%d]: %lf\t",i, j, multRA[i][j]);
	// 	}
	// 	printf("\n");
	// }
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{	
			res[i] += multRA[i][j] * dW[j];
		}
	}
}
