subroutine timestep_caps(s, result)
    type(star_info), pointer :: s
    integer, intent(inout) :: result
    real(dp) :: mass_ratio, dt_at_dMBH_cap, mass_ratio_caps_on, dt_cap_early, dMBH_div_Menv_max, dMBH_div_MBH_max, floor_factor
    logical :: caps_on
    mass_ratio = M_BH / s% mstar
    mass_ratio_caps_on = s% x_ctrl(11)
    dt_cap_early = s% x_ctrl(14)
    dMBH_div_Menv_max = s% x_ctrl(9)
    dMBH_div_MBH_max = s% x_ctrl(10)
    floor_factor = s% x_ctrl(12)

    if (.not. s% lxtra(1) .and. mass_ratio > mass_ratio_caps_on) then
        s% lxtra(1) = .true.
        write(*,'(a,i0,a,f7.4)') '  dM_BH caps on from model ', s% model_number, '  M_BH/Mstar = ', mass_ratio
    end if
    caps_on = s% lxtra(1)

    if (.not. caps_on .and. dt_cap_early > 0d0) then
        s% max_timestep = min(max_timestep_inlist, dt_cap_early)
    else if (caps_on .and. dMBH_div_Menv_max <= 0d0) then
        s% max_timestep = max_timestep_inlist
    end if

    if (caps_on .and. dMBH_div_Menv_max > 0d0) then
        dt_at_dMBH_cap = dMBH_div_Menv_max * s% xmstar * eps * clight**2 / ((1d0 - eps) * s% L_center)
        if (dMBH_div_MBH_max > 0d0) dt_at_dMBH_cap = min(dt_at_dMBH_cap, dMBH_div_MBH_max * M_BH * eps * clight**2 / ((1d0 - eps) * s% L_center))
        s% max_timestep = min(max_timestep_inlist, dt_at_dMBH_cap)
        if (mod(s% model_number, 25) == 0) write(*,'(a,es10.3e2,a)') '  dt-cap: max_timestep = ', s% max_timestep, ' s'
        if (floor_factor > 0d0 .and. s% dt < floor_factor * dt_at_dMBH_cap) then
            write(*,'(a,es9.2e2,a,es9.2e2,a)') '  dt floor: dt = ', s% dt, ' s < ', floor_factor * dt_at_dMBH_cap, ' s: terminating'
            result = terminate
        end if
    end if

    if (s% use_other_wind .and. s% mstar_dot < 0d0) s% max_timestep = min(s% max_timestep, wind_dM_div_M_max * s% mstar / (-s% mstar_dot))
    s% dt_next = min(s% dt_next, s% max_timestep)
end subroutine timestep_caps

subroutine expand_stop(s, result)
    type(star_info), pointer :: s
    integer, intent(inout) :: result
    integer :: n_hist, k
    real(dp) :: dL, dR, L_lim, flat_tol
    flat_tol = s% x_ctrl(6)
    n_hist = max(2, min(expand_hist_max, nint(s% x_ctrl(7))))
    L_lim = s% x_ctrl(8)
    expand_n = expand_n + 1
    expand_L_hist(mod(expand_n - 1, n_hist) + 1) = s% L_center
    expand_R_hist(mod(expand_n - 1, n_hist) + 1) = s% r(1)
    if (expand_n >= n_hist + 50 .and. s% L_center < L_lim * L_Edd_step) then
        dL = 0d0; dR = 0d0
        do k = 1, n_hist
            dL = max(dL, abs(expand_L_hist(k) / s% L_center - 1d0))
            dR = max(dR, abs(expand_R_hist(k) / s% r(1) - 1d0))
        end do
        if (mod(s% model_number, 25) == 0) write(*,'(a,i0,2(a,es9.2e2))') '  expand convergence: model ', s% model_number, '  max|dL/L| = ', dL, '  max|dR/R| = ', dR
        if (dL < flat_tol .and. dR < flat_tol) then
            write(*,'(a,i0,a,i0,a)') '  EXPAND CONVERGED at model ', s% model_number, ': L_center and R flat to x_ctrl(6) over the last ', n_hist, ' steps'
            result = terminate
        end if
    end if
end subroutine expand_stop
