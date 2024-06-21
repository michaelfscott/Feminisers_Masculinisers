#read SLURM_ARRAY_TASK_ID in
args = commandArgs(trailingOnly=TRUE)
array_id=as.numeric(args[1])

source("functions.R")

#common parameters
startFreqF2=0.05
startFreqM2=0.05
end=0.00001

#how many parameter combinations to be run in each array job
per_task=400
#get the start and end IDs for the parameter combinations to be run
START_NUM=((array_id - 1) * per_task) + 1 
END_NUM=( array_id * per_task )

#output table name
output=paste0("../results/simulations/final_frequencies." , START_NUM , "." , END_NUM , ".tsv")

#get the parameter combinations to be used
grid<-read_csv("../results/simulations/parameter_grid.csv") %>%
        slice(START_NUM:END_NUM)

#initialise results
res=list()

#loop over each parameter combination in the grid
for(i in seq(1,nrow(grid))){

#run simulation
tbl<-runSim_Fem_Masc(startFreqF2, startFreqM2, end, Rec=grid$Rec[i], k=grid$k[i], K=grid$K[i], theta=grid$theta[i], delta=grid$delta[i], F2dom=grid$F2dom[i], M2dom=grid$M2dom[i], F2epidom=grid$F2epidom[i], M2epidom=grid$M2epidom[i], min_generations=100, max_generations=10000)

#append only the final frequency to the results table
res[[i]]=tbl[nrow(tbl),]
}

#bind all the final frequencies together
out=do.call(rbind, res)

#write output
write_tsv(out, path=output)

#end of R script


