#!/bin/bash
#SBATCH -A mat751
#SBATCH -t 23:00:00
#SBATCH -o %x-%j.out
##SBATCH -p batch
#SBATCH -p extended
#SBATCH -n 1
#SBATCH --job-name=peo
#SBATCH --nodes=1
#SBATCH --gpus=4
#SBATCH --gpus-per-node=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --threads-per-core=1

module reset
module load PrgEnv-cray/8.4.0
module load amd-mixed/5.7.0
module load craype-accel-amd-gfx90a
module load cce/16.0.1
module load craype/2.7.23
module load cray-fftw
module load rocm/5.7.1

# Path to NAMD3 executable
namd3="/ccs/home/kimsn/program/NAMD_3.0_Source/Linux-x86_64-clang++.mpi/namd3"

# Prefixes for input and output files
#equi_prefix="step4.%d_equilibration"
equi_prefix="step4_equilibration"
prod_prefix="step5_production"
prod_step="step5"

## Uncomment this section to run equilibration steps
#srun "${namd3}" +ppn 28 +ignoresharing  "${equi_prefix}.inp" > "${equi_prefix}.out"



# Running production for 50 steps
for cnt in {5..200};
do
    if [ "$cnt" -eq 1 ]; then
        outputname="${prod_step}_${cnt}"
        sed "s/${prod_prefix}/${outputname}/" "${prod_prefix}.inp" > "${prod_step}_run.inp"
    else
        cntprev=$((cnt - 1))
        equi_last=$(printf "${equi_prefix}" "6")
        inputname="${prod_step}_${cntprev}"
        outputname="${prod_step}_${cnt}"
        sed "s/${equi_last}/${inputname}/" "${prod_prefix}.inp" | \
            sed "s/${prod_prefix}/${outputname}/" > "${prod_step}_run.inp"
    fi
    srun "${namd3}" +ppn 28 +ignoresharing  "${prod_step}_run.inp" > "${outputname}.out"

    # Run the simulation
    #srun --ntasks=31 "${namd3}" +setcpuaffinity +p 8 +devices 0,1 +ignoresharing +idlepoll "${prod_step}_run.inp" > "${outputname}.out"
#    srun "${namd3}" +ppn 7 +ignoresharing  "${prod_step}_run.inp" > "${outputname}_6-7.out"
#    srun "${namd3}" +ppn 14 +ignoresharing  "${prod_step}_run.inp" > "${outputname}_6-14.out"
#    srun "${namd3}" +ppn 21 +ignoresharing  "${prod_step}_run.inp" > "${outputname}_6-21.out"
#    srun "${namd3}" +ppn 28 +ignoresharing  "${prod_step}_run.inp" > "${outputname}_4-28.out"
#    srun "${namd3}" +ppn 35 +ignoresharing  "${prod_step}_run.inp" > "${outputname}_6-35.out"
#    srun "${namd3}" +ppn 42 +ignoresharing  "${prod_step}_run.inp" > "${outputname}_6-42.out"
#    srun "${namd3}" +ppn 49 +ignoresharing  "${prod_step}_run.inp" > "${outputname}_6-49.out"
#    srun "${namd3}" +ppn 56 +ignoresharing  "${prod_step}_run.inp" > "${outputname}_6-56.out"


#    srun "${namd3}" +ppn 7 +ignoresharing  "${prod_step}_run.inp" > "${outputname}_8.out"
#    srun "${namd3}" +ppn 15 +ignoresharing  "${prod_step}_run.inp" > "${outputname}_16.out"
#    srun "${namd3}" +ppn 32 +ignoresharing  "${prod_step}_run.inp" > "${outputname}_32-2.out"
#    srun "${namd3}" +ppn 47 +ignoresharing  "${prod_step}_run.inp" > "${outputname}_48.out"
#    srun "${namd3}" +ppn 55 +ignoresharing  "${prod_step}_run.inp" > "${outputname}_6-56.out"
#    srun "${namd3}" +ppn 111 +ignoresharing  "${prod_step}_run.inp" > "${outputname}_112.out"
done

