#!/bin/bash

#Because there are thousands of simulations, we run them in parallel using an array job. 
#We must first create a grid of parameters to use.
#Then parts of this parameter grid are submitted in parallel in an array job
#Finally the results are combined to give a single table of outputs. 
#The following submits the relevant jobs for the pre-amble (create of parameter grid), array job, and cleanup (combining tables), in the correct order with dependencies on the previous job finishing. 

create_parameter_grid_id=$(sbatch --parsable create_fig3_data1.sub)

run_simulations_array_id=$(sbatch --parsable --dependence=afterok:$create_parameter_grid_id create_fig3_data2.sub)

combine_tables_id=$(sbatch --parsable --dependency=afterok:$run_simulations_array_id create_fig3_data3.sub)

