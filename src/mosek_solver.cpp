//
//  mosek_solver.cpp
//  hmm_aw
//
//  Created by Yukun Chen on 3/8/17.
//  Copyright © 2017 cyk. All rights reserved.
//
#include "mosek.h"
#include "utils/blas_utils.h"
#include "mosek_solver.h"

static MSKenv_t env = NULL;

#include <utility>
#include <map>
using std::pair;
using std::make_pair;
using std::map;
static map< pair<int, int>, MSKtask_t > task_mapper;

void solver_setup() {
  MSKrescodee r;
  /* Create the mosek environment. */
  r = MSK_makeenv(&env, NULL);
  r = MSK_checkoutlicense(env, MSK_FEATURE_PTS);
  if ( r != MSK_RES_OK ) {
    /* In case of an error print error code and description. */
    char symname[MSK_MAX_STR_LEN];
    char desc[MSK_MAX_STR_LEN];
    
    printf("An error occurred during setup.\n");
    MSK_getcodedesc (r,
                     symname,
                     desc);
    printf("Error %s - '%s'\n",symname,desc);
    exit(-1);
  }
}

void solver_release() {
  /*
   size_t i;
   for (i=0; i<task_seq_size; ++i)
   if (task_seq[i] != NULL) MSK_deletetask(&task_seq [i]);
   */
  for (map< pair<int, int>, MSKtask_t >::iterator it=task_mapper.begin(); it!=task_mapper.end(); ++it) MSK_deletetask(&(it->second));
  MSK_deleteenv(&env);
}

double match_by_distmat(int n, int m, double *C, double *wX, double *wY,\
                           /** OUT **/ double *x, /** OUT **/ double *lambda) {
  
  const MSKint32t numvar = n * m,
  numcon = n + m - 1;
  MSKtask_t    *p_task;
  MSKrescodee r = MSK_RES_OK;
  MSKint32t    i,j;
  double fval = 0.0;

  
  if (task_mapper.find(make_pair (n, m)) == task_mapper.end()) {
    task_mapper[make_pair (n, m)] = NULL;
  }
  p_task = &task_mapper[make_pair (n, m)];
  
  if (*p_task == NULL) {
    MSKint32t *asub;
    double ones[2] = {1.0, 1.0};
    asub = (MSKint32t *) malloc(2*m*n* sizeof(MSKint32t));
    for (j=0; j<m; ++j) {
      for (i=0; i<n; ++i) {
        asub[2*(i +j*n)] = i;
        asub[2*(i +j*n) +1] = n+j;
      }
    }
    
    /* Create the optimization task. */
    r = MSK_maketask(env,numcon,numvar,p_task);
    if ( r==MSK_RES_OK ) {
      //r = MSK_linkfunctotaskstream(*p_task,MSK_STREAM_LOG,NULL,printstr);
    }
    
    r = MSK_appendcons(*p_task,numcon);
    r = MSK_appendvars(*p_task,numvar);
    
    for (j=0; j<numvar && r == MSK_RES_OK; ++j) {
      r = MSK_putvarbound(*p_task,
                          j,           /* Index of variable.*/
                          MSK_BK_LO,      /* Bound key.*/
                          0.0,      /* Numerical value of lower bound.*/
                          +MSK_INFINITY);     /* Numerical value of upper bound.*/
      
      i = (j >= (m-1)*n) ? 1 : 2;
      r = MSK_putacol(*p_task,
                      j,           /* Index of variable.*/
                      i,           /* Number of non-zeros in column j.*/
                      asub+j*2,
                      ones);
    }
    free(asub);
    
    if (r == MSK_RES_OK) {
      r = MSK_putobjsense(*p_task, MSK_OBJECTIVE_SENSE_MINIMIZE);
    }
    /* Disable presolve: may lead to minor improvement */
    // r = MSK_putintparam(task, MSK_IPAR_PRESOLVE_USE, MSK_PRESOLVE_MODE_OFF);
    /* set network flow problem */
    r = MSK_putintparam(*p_task, MSK_IPAR_OPTIMIZER,  MSK_OPTIMIZER_NETWORK_PRIMAL_SIMPLEX);
    /* disable multi-threads */
    r = MSK_putintparam(*p_task, MSK_IPAR_NUM_THREADS, 1);
  } else {
  }
  
  /* modify an existing task and re-optimize */
  for (j=0; j<numvar && r == MSK_RES_OK; ++j) {
    r = MSK_putcj(*p_task,j,C[j]);
  }
  
  for (i=0; i<n && r==MSK_RES_OK; ++i)
    r = MSK_putconbound(*p_task,
                        i,
                        MSK_BK_FX,
                        wX[i],
                        wX[i]);
  for (i=0; i<m-1 && r==MSK_RES_OK; ++i)
    r = MSK_putconbound(*p_task,
                        i+n,
                        MSK_BK_FX,
                        wY[i],
                        wY[i]);
  
  
  if ( r==MSK_RES_OK )
  {
    MSKrescodee trmcode;
    
    /* Run optimizer */
    r = MSK_optimizetrm(*p_task,&trmcode);
    
    /* Print a summary containing information
     about the solution for debugging purposes. */
    //MSK_solutionsummary (task,MSK_STREAM_LOG);
    if ( r==MSK_RES_OK ) {
      MSKsolstae solsta;
      r = MSK_getsolsta (*p_task,
                         MSK_SOL_BAS,
                         &solsta);
      switch(solsta)
      {
        case MSK_SOL_STA_OPTIMAL:
        case MSK_SOL_STA_NEAR_OPTIMAL:
        {
          //double *xx = (double*) calloc(numvar,sizeof(double));
          MSK_getprimalobj(*p_task, MSK_SOL_BAS, &fval);
          if ( x )
          {
            MSK_getxx(*p_task,
                      MSK_SOL_BAS,    /* Request the basic solution. */
                      x);
            //printf("Optimal primal solution\n");
            //for(j=0; j<numvar; ++j) printf("x[%d]: %e\n",j,xx[j]);
            
            //free(x);
          }
          
          if (lambda)
          {
            MSK_gety (*p_task,
                      MSK_SOL_BAS,    /* Request the dual solution: be careful about exact +- of variables */
                      lambda);

          }
          else
            r = MSK_RES_ERR_SPACE;
          
          break;
        }
        case MSK_SOL_STA_DUAL_INFEAS_CER:
        case MSK_SOL_STA_PRIM_INFEAS_CER:
        case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
        case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
          printf("Primal or dual infeasibility certificate found.\n");
          break;
        case MSK_SOL_STA_UNKNOWN:
        {
          char symname[MSK_MAX_STR_LEN];
          char desc[MSK_MAX_STR_LEN];
          
          /* If the solutions status is unknown, print the termination code
           indicating why the optimizer terminated prematurely. */
          
          MSK_getcodedesc(trmcode,
                          symname,
                          desc);
          
          printf("The solution status is unknown.\n");
          printf("The optimizer terminitated with code: %s\n",symname);
          break;
        }
        default:
          printf("Other solution status.\n");
          break;
      }
    }
  }
  
  //  MSK_deletetask(&task);
  
  return fval;
}

