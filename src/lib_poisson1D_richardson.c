/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
//all v.p.
void eig_poisson1D(double* eigval, int *la){
    for(int j = 1; j<=*la ; ++j){
         eigval[j-1] = 2*(1-cos((M_PI*j)/ (*la+1)));
     }
}
//vp max
double eigmax_poisson1D(int *la){
  return 2*(1-cos((M_PI*(*la))/ (*la+1)));

}
//vp min
double eigmin_poisson1D(int *la){
    return 2*(1-cos((M_PI* 1)/ (*la+1)));
}
//best alpha
double richardson_alpha_opt(int *la){
    return 2/(eigmin_poisson1D(la)+eigmax_poisson1D(la));
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
    
    double *AX = (double *) calloc((*la), sizeof(double));
    //AX
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1, AB, *lab, X, 1, 0, AX, 1);
    // erreur relatif ||r||/||b||
    double er_rel = relative_forward_error(AX,RHS,la);
 
    //Richardson
    while (er_rel > *tol && *nbite < *maxit) {
        
        for (int i = 0; i < *la; ++i) {
            X[i] = X[i] + (*alpha_rich) * (RHS[i]-AX[i]);
        }
        
        er_rel = relative_forward_error(AX,RHS,la);
        resvec[*nbite] = er_rel;
        
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1, AB, *lab, X, 1, 0, AX, 1);
        
        *nbite += 1;
    }

    free(AX);

}
//D^-1
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
    
    for (int i=0;i<(*la);++i){
        MB[i*(*lab)+ *kv]=1/AB[i*(*lab)+ *kv];
    }
    
}

void extract_MB_gauss_seidel_tridiag(double *AB, double **MB, int *lab, int *la,int *ku, int*kl, int *kv){
    //*lab=3,*la=n
    int num=(*la)*(*la);
    *MB = (double *) calloc((num), sizeof(double));
    int i,j,a,x;
    int dif;
    for (i=0;i<(*la);++i){
        (*MB)[(*la)*i]=1/AB[i*(*lab)+ *kv];
        
    }
   
    for (i=0;i<(*la);++i){
        x=(*la)-i;
        a=i;
        for(j=1;j<x;++j){
            
            (*MB)[(*la)*i+j]=(*MB)[(*la)*i+j-1]*(-1)* AB[a*(*lab)+2]/AB[(a+1)*(*lab)+1];
            a+=1;
           
        }
       
    }
    (*lab)=(*la);

}
    

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
    int lab_AB=3;
    int kl_AB=1;
    int ku_AB=1;
    double *AX = (double *) calloc((*la), sizeof(double));
    double *r = (double *) calloc((*la), sizeof(double));
    //AX
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, kl_AB, ku_AB, 1, AB, lab_AB, X, 1, 0, AX, 1);
   
    // erreur relatif ||r||/||b||
    double er_rel = relative_forward_error(AX,RHS,la);
 
    //Richardson
    while (er_rel > *tol && *nbite < *maxit) {
        //b-AX
        for(int i=0;i<*la;++i){
            r[i]=RHS[i]-AX[i];        }
        
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1, MB, *lab, r, 1, 1, X, 1);
         
        
        er_rel = relative_forward_error(AX,RHS,la);
        resvec[*nbite] = er_rel;
        
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, kl_AB, ku_AB, 1, AB, lab_AB, X, 1, 0, AX, 1);
        
        *nbite += 1;
    }

    free(AX);
    free(r);

}

