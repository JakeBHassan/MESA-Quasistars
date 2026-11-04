real(dp), save :: max_timestep_inlist = 0d0
real(dp), parameter :: wind_dM_div_M_max = 1d-3

integer, parameter :: expand_hist_max = 200
real(dp), save :: expand_L_hist(expand_hist_max) = 0d0, expand_R_hist(expand_hist_max) = 0d0
integer, save :: expand_n = 0
