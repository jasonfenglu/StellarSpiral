//$Id: gsl_wrapper.c 122 2011-07-22 09:30:57Z ccfeng $
#include <math.h>
#include <stdio.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_ellint.h>

 #ifdef __cplusplus
    extern"C" {
 #endif
void wrapper_laplace_coei(double *m,double *beta,double *out);
void wrapper_dlaplace_coei(double *m,double *beta,double *out);

void wrapper_ellipticK(double *K,double *out);
void wrapper_ellipticE(double *K,double *out);

 #ifdef __cplusplus
 }
 #endif

//
// Pass data to FORTRAN
//

void wrapper_ellipticK(double *K,double *out){
	*out = gsl_sf_ellint_Kcomp(pow(*K,0.5),'GSL_PRECPDOUBLE');
}

void wrapper_ellipticE(double *K,double *out){
	*out = gsl_sf_ellint_Ecomp(pow(*K,0.5),'GSL_PRECPDOUBLE');
}


//
// Here is for laplace coefficient integral
//

typedef struct laplace_para{
	double s;		//These three define general form
	double n;
	double h;
	double beta;		//this is the input
} Laplace_para;

double laplace_coefi_general_befor (double x, void * params) {
  Laplace_para *laplace_para = (Laplace_para *) params;
  double s = laplace_para->s;
  double h = laplace_para->h;
  double beta = laplace_para->beta;

  double ans =	2.0/M_PI*(pow(gsl_sf_cos(x)-beta,h))
  /pow(1+pow(beta,2)-2*beta*gsl_sf_cos(x),s);
  return ans;

}

double laplace_coefi(double n,double x){
  double result, error;
  double s = 0.5;
  Laplace_para para = {s, n, 0, x};


  gsl_function F;
  F.function = &laplace_coefi_general_befor;
  F.params = &para;

  //initial QAWO adaptive integration for oscillatory functions
  gsl_integration_qawo_table *qawo_table 
  = gsl_integration_qawo_table_alloc(n, M_PI,GSL_INTEG_COSINE,100);

  //initial gsl integral workspace
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

  //running qawo integral
  gsl_integration_qawo(&F, 0, 1e-13,1e-13,1000,w,qawo_table,&result,&error);

  //deallocate
  gsl_integration_workspace_free (w);
  gsl_integration_qawo_table_free(qawo_table);
  


  //old way to integrate, but can't converge
  /*
  gsl_integration_qags (&F, 0, M_PI, 0, 1e-12, 1000,
                        w, &result, &error);
  */
  return result;
}

double dlaplace_coefi(double n,double x){
  double result, error;
  double s = 1.5;
  Laplace_para para = {s, n, 1, x};

  gsl_function F;
  F.function = &laplace_coefi_general_befor;
  F.params = &para;

  //initial QAWO adaptive integration for oscillatory functions
  gsl_integration_qawo_table *qawo_table 
  = gsl_integration_qawo_table_alloc(n, M_PI,GSL_INTEG_COSINE,100);

  //initial gsl integral workspace
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

  //running qawo integral
  gsl_integration_qawo(&F, 0, 1e-13,1e-13,1000,w,qawo_table,&result,&error);

  //deallocate
  gsl_integration_workspace_free (w);
  gsl_integration_qawo_table_free(qawo_table);


  //old way to integrate, but can't converge
  /*
  gsl_integration_qags (&F, 0, M_PI, 0, 1e-12, 1000,
                        w, &result, &error);
  */

  //result = 0;

  return result;
}

//
// Here is for connecting to FORTRAN
//
void wrapper_laplace_coei(double *m,double *beta,double *out){
	*out = laplace_coefi(*m,*beta);
}

void wrapper_dlaplace_coei(double *m,double *beta,double *out){
	*out = dlaplace_coefi(*m,*beta);
}
