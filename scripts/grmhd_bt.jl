#using ClusterManagers
#using Distributed
#pids1 = addprocs(SlurmManager(40), partition="amdpreq")
#pids2 = addprocs(SlurmManager(40), partition="ehtq")
#pids3 = addprocs(SlurmManager(40), partition="ehtq")
#pids4 = addprocs(SlurmManager(40), partition="ehtq")

include(joinpath(@__DIR__,"parsim.jl"))

pids1 = workers()
parsimorder(pids1, "../data/sims/grmhd_runs/dataset000_timeshifted_3599_LO_LMTcal_JCMTcal_regionIII", "../_research/grmhd_runs_bt/dataset000_timeshifted_3599_LO_LMTcal_JCMTcal_regionIII", 1:4)
parsimorder(pids1, "../data/sims/grmhd_runs/dataset000_3599_LO_LMTcal_JCMTcal_regionIII", "../_research/grmhd_runs_bt/dataset000_3599_LO_LMTcal_JCMTcal_regionIII", 1:4)
parsimorder(pids1, "../data/sims/grmhd_runs/dataset010_3599_LO_LMTcal_JCMTcal_regionIII", "../_research/grmhd_runs_bt/dataset010_3599_LO_LMTcal_JCMTcal_regionIII", 1:4)
parsimorder(pids1, "../data/sims/grmhd_runs/dataset010_masked_fluxscl_3599_LO_LMTcal_JCMTcal_regionIII", "../_research/grmhd_runs_bt/dataset010_masked_fluxscl_3599_LO_LMTcal_JCMTcal_regionIII", 1:4)



rmprocs(workers()...)
