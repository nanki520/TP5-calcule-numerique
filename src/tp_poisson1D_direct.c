/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include <time.h>

#define TRF 0
#define TRI 1
#define SV 2

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
    int IMPLEM = 0;
    
    struct timespec t1_tr, t2_tr, t1_sv, t2_sv;
    int maxite= 55;
    for(int x = 10; x < 100000 ; x = x+1000){
        for (int ite = 0 ; ite<maxite ; ite++){
            int ierr;
            int jj;
            int nbpoints, la;
            int ku, kl, kv, lab;
            int *ipiv;
            int info = 1;
            int NRHS;
            double T0, T1;
            double *RHS, *RHS_old, *EX_SOL, *X;
            double **AAB;
            double *AB;
            double *AB_bond;
            double *Y, *Y1;
            
            double relres;
            
            if (argc == 2) {
                IMPLEM = atoi(argv[1]);
            } else if (argc > 2) {
                perror("Application takes at most one argument");
                exit(1);
            }
            
            NRHS=1;
            nbpoints=x;
            la=nbpoints-2;
            T0=-5.0;
            T1=5.0;
            
            printf("--------- Poisson 1D ---------\n\n");
            RHS=(double *) malloc(sizeof(double)*la);
            RHS_old=(double *) malloc(sizeof(double)*la);
            EX_SOL=(double *) malloc(sizeof(double)*la);
            X=(double *) malloc(sizeof(double)*la);
            
            // TODO : you have to implement those functions
            set_grid_points_1D(X, &la);
            set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
            set_dense_RHS_DBC_1D(RHS_old,&la,&T0,&T1);
            set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
            
            //write_vec(RHS, &la, "RHS.dat");
            //write_vec(EX_SOL, &la, "EX_SOL.dat");
            //write_vec(X, &la, "X_grid.dat");
            
            kv=1;
            ku=1;
            kl=1;
            lab=kv+kl+ku+1;
            
            AB = (double *) malloc(sizeof(double)*lab*la);
            AB_bond= (double *) malloc(sizeof(double)*lab*la);
            
            set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
            set_GB_operator_colMajor_poisson1D(AB_bond, &lab, &la, &kv);
            
            //write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");
            
            //printf("Solution with LAPACK\n");
            info=0;
            ipiv = (int *) calloc(la, sizeof(int));
            Y= (double *) calloc(la, sizeof(double));
            Y1= (double *) calloc(la, sizeof(double));
            
            
            /* LU Factorization */
            if (IMPLEM == TRF) {
                clock_gettime(CLOCK_MONOTONIC_RAW, &t1_tr);
                dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
                /* Solution (Triangular) */
                if (info==0){
                    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info,1);
                    if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
                    clock_gettime(CLOCK_MONOTONIC_RAW, &t2_tr);
                }else{
                    printf("\n INFO = %d\n",info);
                }
            }
            
            /* LU for tridiagonal matrix  (can replace dgbtrf_) */
            if (IMPLEM == TRI) {
                dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
                /* Solution (Triangular) */
                if (info==0){
                    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info,1);
                    if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
                    clock_gettime(CLOCK_MONOTONIC_RAW, &t2_tr);
                }else{
                    printf("\n INFO = %d\n",info);
                }
            }
            
            
            /* It can also be solved with dgbsv */
            if (IMPLEM == SV) {
                
                clock_gettime(CLOCK_MONOTONIC_RAW, &t1_sv);
                // TODO : use dgbsv
                dgbsv_( &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS,&la, &info);
                clock_gettime(CLOCK_MONOTONIC_RAW, &t2_sv);
            }
            
            //write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
            //write_xy(RHS, X, &la, "SOL.dat");
            
            /* Relative forward error */
            relres = relative_forward_error(RHS, EX_SOL,&la);
            printf("\nThe relative forward error is relres = %e\n",relres);
            
            /*
             Validation cblas_dgbmv
             
            cblas_dgbmv(CblasColMajor,CblasNoTrans,la,la,kl,ku,1,AB_bond,lab,RHS,1,1,Y,1); //AX
            cblas_dgbmv(CblasColMajor,CblasNoTrans,la,la,kl,ku,1,AB_bond,lab,EX_SOL,1,1,Y1,1); //b
            double error_relative=relative_forward_error(Y, Y1, &la);
            printf("error_relative = %f\n",error_relative);
            */
            
            free(AB_bond);
            free(Y);
            free(Y1);
            
            free(RHS);
            free(EX_SOL);
            free(X);
            free(AB);
            printf("\n\n--------- End -----------\n");
        }
        
        if (IMPLEM == TRF) {
            FILE *file;
            file = fopen("TRF", "a");
            double elapsed = (double)(t2_tr.tv_nsec - t1_tr.tv_nsec)/maxite;
            fprintf(file, "%d %f\n", x,elapsed);
            fclose(file);
        }
        
        if (IMPLEM == SV) {
            FILE *file;
            file = fopen("SV", "a");
            double elapsed = (double)(t2_sv.tv_nsec - t1_sv.tv_nsec)/maxite;
            fprintf(file, "%d %f\n", x,elapsed);
            fclose(file);
        }
    }
}

