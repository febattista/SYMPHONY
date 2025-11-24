#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include <time.h>

#include "symphony.h"
#include "../include/sym_types.h"
#include "../include/sym_master.h"

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
   // sym_set_int_param(env_warm, "max_active_nodes", 1);

   /* Do a first solve */
   sym_solve(env_warm);
   double dual_bound_from_tree, dual_bound_from_df;
   sym_get_obj_val(env_warm, &dual_bound_from_tree);

   sym_build_dual_func(env_warm);

   /* Change the right-hand side */
   double new_rhs = -2000.0;

   sym_evaluate_dual_function(env_warm, &new_rhs, 1, &dual_bound_from_df);

   sym_get_dual_bound_from_tree(env_warm, &new_rhs, 1, &dual_bound_from_tree);

   // Memory clean-ups
   sym_close_environment(env_warm);

   return 0;
}  
