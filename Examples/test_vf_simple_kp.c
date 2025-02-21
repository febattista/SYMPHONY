#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include <time.h>

#include "symphony.h"
// feb223
// These paths should be changed
// I still don't figure out how to use makefile to 
// include these files automatically
#include "../include/sym_types.h"
#include "../include/sym_master.h"

void generate_rand_array(int m, int l, int u, double *res) {
  for (int i = 0; i < m; i++) {
    res[i] = (rand() % (u - l + 1)) + l;
  }
}

void set_rhs(sym_environment *env, double *rhs, int m){
   for (int i = 0; i < m; i++){
      switch (env->mip->sense[i]){
         case 'E' :
            sym_set_row_upper(env, i, rhs[i]);
            sym_set_row_lower(env, i, rhs[i]);
            break;
         case 'G':
            sym_set_row_lower(env, i, rhs[i]);
            break;
         case 'L':
            sym_set_row_upper(env, i, rhs[i]);
            break;
         default:
            printf("Error\n");
      }
   }
}

double *linspace(double a, double b, int pieces){

   double *res = (double *)malloc(sizeof(double) * pieces);
   double step = fabs(a + b)/((double)pieces);

   res[0] = a;
   res[pieces - 1] = b;

   for (int i = 1; i < pieces - 1; i++){
      res[i] = res[i - 1] + step;
   }

   return res;
}

int main(int argc, char **argv)
{   
   
   int termcode;
   int numTrain = 10;
   int numTests = 10;
   double warmObjVal, coldObjVal, dualFuncObj;

   sym_environment *env_warm = sym_open_environment(); 

   sym_parse_command_line(env_warm, argc, argv); 

   sym_load_problem(env_warm);

   sym_set_int_param(env_warm, "verbosity", 10);

   sym_environment *env_cold = sym_open_environment(); 

   sym_parse_command_line(env_cold, argc, argv); 

   sym_load_problem(env_cold);

   // sym_set_int_param(env_cold, "verbosity", -2);

   // sym_set_int_param(env_warm, "keep_warm_start", TRUE);
   // sym_set_int_param(env_warm, "keep_dual_function_description", TRUE);
   // sym_set_int_param(env_warm, "should_use_rel_br", FALSE);
   // sym_set_int_param(env_warm, "use_hot_starts", FALSE);
   // sym_set_int_param(env_warm, "should_warmstart_node", TRUE);
   // sym_set_int_param(env_warm, "sensitivity_analysis", TRUE);
   // sym_set_int_param(env_warm, "sensitivity_rhs", true);
   // sym_set_int_param(env_warm, "sensitivity_bounds", TRUE);
   // sym_set_int_param(env_warm, "set_obj_upper_lim", FALSE);
   sym_set_int_param(env_warm, "do_primal_heuristic", FALSE);
   sym_set_int_param(env_warm, "prep_level", -1);
   // sym_set_int_param(env_warm, "tighten_root_bounds", FALSE);
   // sym_set_int_param(env_warm, "max_sp_size", 100);
   // sym_set_int_param(env_warm, "do_reduced_cost_fixing", FALSE);
   // sym_set_int_param(env_warm, "generate_cgl_cuts", FALSE);
   // sym_set_int_param(env_warm, "max_active_nodes", 1);
   // sym_set_int_param(env_warm, "max_presolve_iter", 0);
   // sym_set_int_param(env_warm, "limit_strong_branching_time", 0);

   sym_solve(env_warm);
   exit(0);

   //----------------------
   // linspace numpy-like 
   //----------------------
   double a = -56.5, b = 5, zerotol = 1e-7;
   int pieces = 200;
   double *zeta_lst_large = linspace(-12.5, -10, pieces);
   double *zeta_lst_small = linspace(-15, 5, pieces);
   double *rvf_lst  = (double*)malloc(sizeof(double) * pieces);
   double *df_lst   = (double*)malloc(sizeof(double) * pieces);
   double *rvf_lst_small  = (double*)malloc(sizeof(double) * pieces);
   double *df_lst_small   = (double*)malloc(sizeof(double) * pieces);
   double *rhs = (double*)malloc(sizeof(double) * 1);

   int num_objs = 1;

   // printf("===============================\n");
   // printf("      ZETA LST LARGE\n");
   // printf("===============================\n");
   // printf("zeta_lst_large = [");
   // for (int i = 0; i < pieces; i++){
   //    printf("%.10f, ", zeta_lst_large[i]);
   // }
   // printf("]\n");

   // printf("===============================\n");
   // printf("      ZETA LST SMALL\n");
   // printf("===============================\n");
   // printf("zeta_lst_small = [");
   // for (int i = 0; i < pieces; i++){
   //    printf("%.10f, ", zeta_lst_small[i]);
   // }
   // printf("]\n");

   // // compute rvf
   // printf("===============================\n");
   // printf("          RVF\n");
   // printf("===============================\n");
   printf("rvf_lst = [");
   for (int i = 0; i < pieces; i++){
      rhs = zeta_lst_large + i;
      printf("RHS: %.5f\n", *rhs);
      set_rhs(env_cold, rhs, num_objs);

      sym_solve(env_cold);
      sym_get_obj_val(env_cold, rvf_lst + i);
      printf("%.10f, ", rvf_lst[i]);
   }
   printf("]\n");

   // // compute rvf
   // printf("===============================\n");
   // printf("          RVF\n");
   // printf("===============================\n");
   // printf("rvf_lst_small = [");
   // for (int i = 0; i < pieces; i++){
   //    rhs = zeta_lst_small + i;
   //    set_rhs(env_cold, rhs, num_objs);

   //    sym_solve(env_cold);
   //    sym_get_obj_val(env_cold, rvf_lst_small + i);
   //    printf("%.10f, ", rvf_lst_small[i]);
   // }
   // printf("]\n");
   
   // Build and evaluate DF
   *rhs = 40.0/9.0;
   printf("===============================\n");
   printf("          DF: %.3f\n", *rhs);
   printf("===============================\n");
   set_rhs(env_warm, rhs, num_objs);
   sym_warm_solve(env_warm);
   sym_build_dual_func(env_warm);
   *rhs = -55.5;
   printf("===============================\n");
   printf("          DF: %.3f\n", *rhs);
   printf("===============================\n");
   set_rhs(env_warm, rhs, num_objs);
   sym_warm_solve(env_warm);
   sym_build_dual_func(env_warm);
   *rhs = -73.0/6.0;
   printf("===============================\n");
   printf("          DF: %.3f\n", *rhs);
   printf("===============================\n");
   set_rhs(env_warm, rhs, num_objs);
   sym_warm_solve(env_warm);
   sym_build_dual_func(env_warm);
   *rhs = -11;
   printf("===============================\n");
   printf("          DF: %.3f\n", *rhs);
   printf("===============================\n");
   set_rhs(env_warm, rhs, num_objs);
   sym_warm_solve(env_warm);
   sym_build_dual_func(env_warm);
   printf("df_lst = [");
   for (int i = 0; i < pieces; i++){
      rhs = zeta_lst_large + i;
      printf("RHS: %.5f\n", *rhs);
      sym_evaluate_dual_function(env_warm, rhs, num_objs, df_lst + i);
      printf("%.10f, ", df_lst[i]);
   }
   printf("]\n");
   // printf("===============================\n");
   // printf("df_lst_small = [");
   // for (int i = 0; i < pieces; i++){
   //    rhs = zeta_lst_small + i;
   //    sym_evaluate_dual_function(env_warm, rhs, num_objs, df_lst_small + i);
   //    printf("%.10f, ", df_lst_small[i]);
   // }
   // printf("]\n");

   // Memory clean-ups
   sym_close_environment(env_warm);
   sym_close_environment(env_cold);
   // free(zetas);

   return 0;
}  

