integer, parameter :: first_solve = 0, trial_pending = 1, trial_applied = 2, step_accepted = 3
integer,  save :: iteration_phase = first_solve
real(dp), save :: L_trial = 0d0

real(dp), parameter :: r_lo_div_R_i = 1.01d0, r_hi_div_R_i = 1.06d0
integer,  parameter :: max_iterations = 15
real(dp), parameter :: tolerance = 1d-2
real(dp), parameter :: first_trial_dlnL = 2d-3
real(dp), parameter :: trial_dlnL_max = 0.25d0
real(dp), parameter :: start_dlnL_max = 0.3d0

integer,  save :: n_iterations = 0
real(dp), save :: lnL_trial_prev = 0d0, f_trial_prev = 0d0
logical,  save :: have_prev_trial = .false.

real(dp), save :: L_accepted = 0d0
logical,  save :: have_L_accepted = .false.
real(dp), save :: L_accepted_prev = 0d0, dlnL_dt = 0d0
real(dp), save :: L_accepted_prev2 = 0d0, dlnL_dt_prev = 0d0
logical,  save :: have_dlnL_dt = .false.
real(dp), save :: slope_prev_step = 0d0, dt_prev_step = 0d0
logical,  save :: have_slope = .false.

real(dp), save :: L_sat_i_last = 0d0, L_conv_i_last = 0d0
