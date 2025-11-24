#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include <time.h>

#include "symphony.h"
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

   sym_environment *env_warm = sym_open_environment(); 

   sym_parse_command_line(env_warm, argc, argv); 

   sym_load_problem(env_warm);

   sym_set_int_param(env_warm, "verbosity", -2);

   // feb223
   // Those are the parameters to be set in order to 
   // keep the branch-and-bound tree valid for RHS changes
   // (E.g. Cuts and reduced cost fixing are not RHS-invariant)
   sym_set_int_param(env_warm, "keep_warm_start", TRUE);
   sym_set_int_param(env_warm, "keep_dual_function_description", TRUE);
   sym_set_int_param(env_warm, "should_use_rel_br", FALSE);
   sym_set_int_param(env_warm, "use_hot_starts", FALSE);
   sym_set_int_param(env_warm, "should_warmstart_node", TRUE);
   sym_set_int_param(env_warm, "sensitivity_analysis", TRUE);
   sym_set_int_param(env_warm, "sensitivity_rhs", true);
   sym_set_int_param(env_warm, "sensitivity_bounds", TRUE);
   sym_set_int_param(env_warm, "set_obj_upper_lim", FALSE);
   sym_set_int_param(env_warm, "do_primal_heuristic", FALSE);
   sym_set_int_param(env_warm, "prep_level", -1);
   sym_set_int_param(env_warm, "tighten_root_bounds", FALSE);
   sym_set_int_param(env_warm, "max_sp_size", 100);
   sym_set_int_param(env_warm, "do_reduced_cost_fixing", FALSE);
   sym_set_int_param(env_warm, "generate_cgl_cuts", FALSE);
   sym_set_int_param(env_warm, "max_active_nodes", 1);

   int i = 0;
   int num_samples = 10;
   int num_objs = 2;
   int termcode;
   double dualFuncObj = 0;
   double warmObjFunc = 0;
   double infinity = 1e10;

   double **rhss = (double **)malloc(sizeof(double*) * num_objs);
   double *objVals = (double *)malloc(sizeof(double) * num_samples);
   for (i = 0; i < num_objs; i++){
      rhss[i] = linspace(0, 3000, num_samples);
   }

   /************************************************************************\
    * Assert the dual function is strong at the tested RHSs
   \************************************************************************/
   for (i = 0; i < num_samples; i++){
      double *rhs = (double *)malloc(sizeof(double) * num_objs);
      rhs[0] = rhss[0][i];
      rhs[1] = rhss[1][i];

      set_rhs(env_warm, rhs, num_objs);
      
      
      if ((termcode = sym_warm_solve(env_warm)) < 0){
         printf("WARM: PROBLEM INFEASIBLE!\n");
      }

      if (sym_is_proven_optimal(env_warm)){
         
      } else if (sym_is_proven_primal_infeasible(env_warm)){
   
      }

      sym_get_obj_val(env_warm, &warmObjFunc);
      sym_build_dual_func(env_warm);
      sym_evaluate_dual_function(env_warm, rhs, num_objs, &dualFuncObj);
      print_dual_function(env_warm);
      if (sym_is_proven_optimal(env_warm)){
         printf("RHS: [%.2f, %.2f]  -->  PROVEN OPTIMAL!\n", 
            rhs[0], rhs[1]);
         printf(" DF: %.10f  RVF: %.10f\n", dualFuncObj, warmObjFunc);
         objVals[i] = warmObjFunc;
         assert(fabs(dualFuncObj - warmObjFunc) < 1e-5); 
      } else if (sym_is_proven_primal_infeasible(env_warm)){
         printf("RHS: [%.2f, %.2f]  -->  PROVEN INFEASIBLE!\n", rhs[0], rhs[1]);
         printf(" DF: %.10f\n", dualFuncObj, warmObjFunc);
         objVals[i] = infinity;
         assert(dualFuncObj > infinity);
      }

      // print_dual_function(env_warm);
      free(rhs);
   }

   print_dual_function(env_warm);

   /************************************************************************\
    * Assert the dual function remains strong at tested RHSs
   \************************************************************************/
   for (i = 0; i < num_samples; i++){
      double *rhs = (double *)malloc(sizeof(double) * num_objs);
      rhs[0] = rhss[0][i];
      rhs[1] = rhss[1][i];

      sym_evaluate_dual_function(env_warm, rhs, num_objs, &dualFuncObj);
      if (objVals[i] > infinity/10){
         assert(dualFuncObj > infinity);
         printf("RHS: [%.2f, %.2f]  -->  PROVEN INFEASIBLE!\n", rhs[0], rhs[1]);
         printf(" DF: %.10f\n", dualFuncObj);
      } else {
         printf("RHS: [%.2f, %.2f]  -->  PROVEN OPTIMAL!\n", 
            rhs[0], rhs[1]);
         printf(" DF: %.10f  RVF: %.10f\n", dualFuncObj, objVals[i]);
         assert(fabs(dualFuncObj - objVals[i]) < 1e-5); 
      }
     
      free(rhs);
   }

   // Memory clean-ups
   sym_close_environment(env_warm);
   free(objVals);
   for (i = 0; i < num_objs; i++){
      free(rhss[i]);
   }
   free(rhss);

   return 0;
}  
