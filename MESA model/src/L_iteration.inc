real(dp) function starting_L_center(s, dt)
    type(star_info), pointer :: s
    real(dp), intent(in) :: dt
    if (have_L_accepted) then
        starting_L_center = L_accepted
        if (s% use_other_wind .and. have_dlnL_dt .and. dlnL_dt * dlnL_dt_prev < 0d0 .and. L_accepted_prev2 > 0d0) then
            starting_L_center = 0.5d0 * (L_accepted + L_accepted_prev2)
        else if (have_dlnL_dt) then
            starting_L_center = starting_L_center * exp(max(-start_dlnL_max, min(start_dlnL_max, dlnL_dt * dt)))
        end if
    else if (s% L_center > 0d0) then
        starting_L_center = s% L_center
    else
        starting_L_center = L_Edd_step
    end if
end function starting_L_center

subroutine zones_near_R_i(s, k_out, k_in)
    type(star_info), pointer :: s
    integer, intent(out) :: k_out, k_in
    integer :: k
    k_in = s% nz - 2
    do while (s% r(k_in) < r_lo_div_R_i * s% R_center)
        k_in = k_in - 1
    end do
    k_out = k_in
    do k = k_in - 1, 1, -1
        if (s% r(k) > r_hi_div_R_i * s% R_center) exit
        k_out = k
    end do
end subroutine zones_near_R_i

subroutine check_L_center(s, result)
    type(star_info), pointer :: s
    integer, intent(inout) :: result
    real(dp) :: L_sat_i, L_conv_i, f, lnL, dlnL, slope
    logical :: converged

    result = keep_going
    if (iteration_phase == step_accepted) return

    call base_convection(s, L_sat_i, L_conv_i)
    L_sat_i_last = L_sat_i
    L_conv_i_last = L_conv_i

    f = (L_sat_i - L_conv_i) / s% L_center
    n_iterations = n_iterations + 1
    converged = abs(f) <= tolerance

    write(*,'(a,i3,a,es11.3e2,a,es11.3e2,a,es10.2e2)') '  L iteration ', n_iterations, ': L/L_Edd=', s% L_center / L_Edd_step, '  L_sat/L_conv=', L_sat_i / L_conv_i, '  f=', f

    lnL = log(s% L_center)
    if (have_prev_trial) slope = (f - f_trial_prev) / (lnL - lnL_trial_prev)

    if (converged .or. n_iterations >= max_iterations) then
        if (converged) then
            write(*,'(a,i3,a,es9.2e2,a)') '  L iteration CONVERGED (', n_iterations, ' iterations, |f| = ', abs(f), ')'
            if (have_prev_trial) then
                slope_prev_step = slope
                dt_prev_step = s% dt
                have_slope = .true.
            end if
        else
            write(*,'(a,i3,a,es10.2e2)') '  L iteration HIT CAP (', n_iterations, ' iterations); accepting with f=', f
        end if
        call accept_L_center(s, converged)
        return
    end if

    if (have_prev_trial) then
        dlnL = -f / slope
    else if (have_slope) then
        dlnL = -f / (slope_prev_step * max(0.1d0, min(10d0, s% dt / dt_prev_step)))
    else
        dlnL = 0.5d0 * log(1d0 + f)
        dlnL = max(-first_trial_dlnL, min(first_trial_dlnL, dlnL))
    end if
    dlnL = max(-trial_dlnL_max, min(trial_dlnL_max, dlnL))

    lnL_trial_prev = lnL
    f_trial_prev = f
    have_prev_trial = .true.

    L_trial = s% L_center * exp(dlnL)
    iteration_phase = trial_pending
    result = redo
end subroutine check_L_center

subroutine accept_L_center(s, update_dlnL_dt)
    type(star_info), pointer :: s
    logical, intent(in) :: update_dlnL_dt
    L_accepted = s% L_center
    have_L_accepted = .true.
    iteration_phase = step_accepted
    L_start = s% L_center
    if (update_dlnL_dt .and. L_accepted_prev > 0d0) then
        dlnL_dt_prev = dlnL_dt
        dlnL_dt = log(s% L_center / L_accepted_prev) / s% dt
        have_dlnL_dt = .true.
    end if
end subroutine accept_L_center

subroutine apply_trial_L_center(s)
    type(star_info), pointer :: s
    if (iteration_phase == trial_applied) then
        have_prev_trial = .false.
        L_trial = s% L_center + 0.5d0 * (L_trial - s% L_center)
        write(*,'(a,es14.5e3)') '  L iteration: the solver failed on the trial; halving it to L_center=', L_trial
    end if
    iteration_phase = trial_applied
    s% L_center = L_trial
    s% need_to_setvars = .true.
end subroutine apply_trial_L_center

subroutine reset_L_iteration()
    iteration_phase = first_solve
    n_iterations = 0
    have_prev_trial = .false.
    L_accepted_prev2 = L_accepted_prev
    L_accepted_prev = L_accepted
end subroutine reset_L_iteration

subroutine apply_boundary_before_solve(id, dt, res)
    integer, intent(in) :: id
    real(dp), intent(in) :: dt
    integer, intent(out) :: res
    type(star_info), pointer :: s
    integer :: ierr
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) then
        res = terminate
        return
    end if
    res = keep_going

    if (core_inserted) then
        call set_inner_boundary(s, dt)
        call star_set_vars(id, dt, ierr)
        if (ierr /= 0) then
            res = terminate
            return
        end if
    end if
    if (iteration_phase == trial_pending .or. iteration_phase == trial_applied) then
        call apply_trial_L_center(s)
    else if (iteration_phase == step_accepted) then
        iteration_phase = first_solve
        n_iterations = 0
        have_prev_trial = .false.
    end if
end subroutine apply_boundary_before_solve
