#!/bin/bash

#mpirun -np 4 -machinefile mpirun_hosts tacc_affinity_07 ../../ks_imp_dyn2/su3_rmd < ./small.in > o07.out &
#mpirun -np 4 -machinefile mpirun_hosts tacc_affinity_16 ../../ks_imp_dyn2/su3_rmd < ./small.in > o16.out &
#mpirun -np 4 -machinefile mpirun_hosts tacc_affinity_25 ../../ks_imp_dyn2/su3_rmd < ./small.in > o25.out &
#mpirun -np 4 -machinefile mpirun_hosts tacc_affinity_34 ../../ks_imp_dyn2/su3_rmd < ./small.in > o34.out &
#wait

#mpirun -np 4 -machinefile mpirun_hosts tacc_affinity_07 ../../ks_imp_dyn2/su3_rmd < ./small.in > o07.out
#mpirun -np 4 -machinefile mpirun_hosts tacc_affinity_16 ../../ks_imp_dyn2/su3_rmd < ./small.in > o16.out
#mpirun -np 4 -machinefile mpirun_hosts tacc_affinity_25 ../../ks_imp_dyn2/su3_rmd < ./small.in > o25.out
#mpirun -np 4 -machinefile mpirun_hosts tacc_affinity_34 ../../ks_imp_dyn2/su3_rmd < ./small.in > o34.out

#mpirun -np 4 -machinefile mpirun_hosts tacc_affinity_07i ../../ks_imp_dyn2/su3_rmd < ./small.in > o07.out &
#mpirun -np 4 -machinefile mpirun_hosts tacc_affinity_16i ../../ks_imp_dyn2/su3_rmd < ./small.in > o16.out &
#mpirun -np 4 -machinefile mpirun_hosts tacc_affinity_25i ../../ks_imp_dyn2/su3_rmd < ./small.in > o25.out &
#mpirun -np 4 -machinefile mpirun_hosts tacc_affinity_34i ../../ks_imp_dyn2/su3_rmd < ./small.in > o34.out &
#wait

mpirun -np 4 -machinefile mpirun_hosts tacc_affinity_07i ../../ks_imp_dyn2/su3_rmd < ./small.in > o07.out
mpirun -np 4 -machinefile mpirun_hosts tacc_affinity_16i ../../ks_imp_dyn2/su3_rmd < ./small.in > o16.out
mpirun -np 4 -machinefile mpirun_hosts tacc_affinity_25i ../../ks_imp_dyn2/su3_rmd < ./small.in > o25.out
mpirun -np 4 -machinefile mpirun_hosts tacc_affinity_34i ../../ks_imp_dyn2/su3_rmd < ./small.in > o34.out

