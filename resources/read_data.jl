using MAT

vars = matread("resources/connectivities_estimates_withsurr_twoAALnets.mat")
DMNlist = vars["DMNlist"]
dmn_results = vars["dmn_results"]
FPNlist = vars["FPNlist"]
fpn_results = vars["fpn_results"]

@show DMNlist
# @show dmn_results
@show FPNlist
# @show fpn_results
