module run_star_extras

use star_lib
use star_def
use const_def
use math_lib
use chem_def

implicit none

real(dp), save :: eps = 0d0
real(dp), save :: beta = 0d0
real(dp), save :: N_rs = 0d0

real(dp), save :: M_BH = 0d0
logical, save :: bh_frozen = .false.

real(dp), save :: L_Edd_step = 0d0
real(dp), save :: Mdot_in = 0d0
real(dp), save :: M_BH_start = 0d0
real(dp), save :: L_start = 0d0
logical,  save :: core_inserted = .false.

include 'L_iteration_def.inc'
include 'timestep_control_def.inc'
include 'other_winds/other_winds_def.inc'

contains

include 'L_iteration.inc'
include 'timestep_control.inc'
include 'other_cgrav.inc'
include 'other_winds/other_winds.inc'
include 'other_winds/other_winds_columns.inc'

subroutine extras_controls(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    real(dp) :: M_BH_initial
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    s% extras_startup => extras_startup
    s% extras_start_step => extras_start_step
    s% extras_check_model => extras_check_model
    s% extras_finish_step => extras_finish_step
    s% how_many_extra_history_columns => how_many_extra_history_columns
    s% data_for_extra_history_columns => data_for_extra_history_columns
    s% how_many_extra_profile_columns => how_many_extra_profile_columns
    s% data_for_extra_profile_columns => data_for_extra_profile_columns
    s% other_cgrav => my_other_cgrav
    s% other_wind => erupt_other_wind

    eps = s% x_ctrl(1)
    beta = s% x_ctrl(3)
    N_rs = s% x_ctrl(2)
    s% job% new_core_mass = s% x_ctrl(4) * s% initial_mass
    s% job% relax_core_years_for_dt = s% x_ctrl(5)
    if (s% x_ctrl(5) <= 0d0) s% job% relax_core_years_for_dt = min(1d0, (s% initial_mass / 1d6)**(-3))
    M_BH_initial = s% job% new_core_mass * Msun
    s% job% core_avg_rho = M_BH_initial / (four_thirds_pi * R_i(M_BH_initial)**3)
    s% job% core_avg_eps = L_Edd(s% initial_mass * Msun) / M_BH_initial
    write(*,'(a,es10.3,a,es10.3,a)') '  BH insertion: new_core_mass = ', s% job% new_core_mass, ' Msun, relax_core_years_for_dt = ', s% job% relax_core_years_for_dt, ' yr (extras_controls)'
    s% use_other_before_struct_burn_mix = .true.
    s% other_before_struct_burn_mix => apply_boundary_before_solve
end subroutine extras_controls

real(dp) function L_sat(s, k)
    type(star_info), pointer :: s
    integer, intent(in) :: k
    L_sat = beta * 4d0 * pi * s% r(k)**2 * s% Peos(k) * s% csound(k)
    L_sat = L_sat * s% gamma1(k) / (4d0 * (s% gamma1(k) - 1d0))
end function L_sat

real(dp) function R_i(M_BH)
    real(dp), intent(in) :: M_BH
    R_i = N_rs * 2d0 * standard_cgrav * M_BH / clight**2
end function R_i

real(dp) function L_Edd(M)
    real(dp), intent(in) :: M
    L_Edd = 4d0 * pi * clight * standard_cgrav * M / 0.34d0
end function L_Edd

subroutine base_convection(s, L_sat_i, L_conv_i)
    type(star_info), pointer :: s
    real(dp), intent(out) :: L_sat_i, L_conv_i
    integer :: k, k_out, k_in
    call zones_near_R_i(s, k_out, k_in)
    L_sat_i = 0d0; L_conv_i = 0d0
    do k = k_out, k_in
        L_sat_i = L_sat_i + L_sat(s, k)
        L_conv_i = L_conv_i + s% L_conv(k)
    end do
    L_sat_i = L_sat_i / dble(k_in - k_out + 1)
    L_conv_i = L_conv_i / dble(k_in - k_out + 1)
end subroutine base_convection

subroutine start_accretion_step(s)
    type (star_info), pointer :: s
    L_Edd_step = L_Edd(s% mstar)
    L_start = starting_L_center(s, s% dt_next)
    Mdot_in = L_start / (eps * clight**2)
    M_BH_start = M_BH
end subroutine start_accretion_step

real(dp) function Mdot_BH()
    Mdot_BH = (1d0 - eps) * Mdot_in
    if (bh_frozen) Mdot_BH = 0d0
end function Mdot_BH

subroutine set_inner_boundary(s, dt_now)
    type(star_info), pointer :: s
    real(dp), intent(in) :: dt_now
    M_BH = M_BH_start + Mdot_BH() * dt_now
    s% xtra(1) = M_BH
    s% M_center = M_BH
    s% xmstar = s% mstar - s% M_center
    s% L_center = L_start
    call do1_relax_R_center(s, R_i(M_BH))
end subroutine set_inner_boundary

subroutine insert_black_hole(id, s, ierr)
    integer, intent(in) :: id
    type (star_info), pointer :: s
    integer, intent(out) :: ierr
    real(dp) :: new_core_mass, core_avg_rho, core_avg_eps
    M_BH = s% M_center
    if (M_BH <= 0d0) M_BH = s% job% new_core_mass * Msun
    call start_accretion_step(s)
    M_BH = M_BH + min(Mdot_BH() * s% dt_next, 0.03d0 * M_BH)
    s% xtra(1) = M_BH
    new_core_mass = M_BH / Msun
    core_avg_rho = M_BH / (four_thirds_pi * R_i(M_BH)**3)
    core_avg_eps = L_start / (new_core_mass * Msun)
    write(*,'(a,es12.4e3,a,f0.3)') '  BH insertion: M_BH = ', new_core_mass, ' Msun,  L_start/L_Edd ', L_start / L_Edd_step
    call star_relax_core(id, new_core_mass, s% job% dlg_core_mass_per_step, s% job% relax_core_years_for_dt, core_avg_rho, core_avg_eps, ierr)
end subroutine insert_black_hole

subroutine do1_relax_R_center(s, new_Rcenter)
    type (star_info), pointer :: s
    real(dp), intent(in) :: new_Rcenter
    real(dp) :: dr3, rp13
    integer :: k
    s% R_center = new_Rcenter
    rp13 = new_Rcenter**3
    do k = s% nz, 1, -1
        dr3 = s% dm(k) / (s% rho(k) * four_thirds_pi)
        s% xh(s% i_lnR, k) = log(rp13 + dr3) * one_third
        rp13 = rp13 + dr3
    end do
end subroutine do1_relax_R_center

subroutine extras_startup(id, restart, ierr)
    integer, intent(in) :: id
    logical, intent(in) :: restart
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    bh_frozen = .not. s% x_logical_ctrl(1)
    if (bh_frozen) write(*,'(a)') '  black hole frozen (x_logical_ctrl(1) = .false.)'
    max_timestep_inlist = s% max_timestep

    if (s% model_number == 0) then
        call insert_black_hole(id, s, ierr)
    else
        M_BH = s% xtra(1)
    end if
end subroutine extras_startup

integer function extras_start_step(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_start_step = keep_going

    call reset_L_iteration()
    call start_accretion_step(s)
    call set_inner_boundary(s, s% dt_next)
    core_inserted = .true.
    write(*,'(a,i0,a,es9.2e2,a,f0.5,a,f0.3)') ' model ', s% model_number, '  dt(s) ', s% dt_next, '  M_BH/Mstar: ', M_BH / s% mstar, '  L_start/L_Edd ', L_start / L_Edd_step
end function extras_start_step

integer function extras_check_model(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    real(dp) :: mass_ratio, mass_ratio_stop
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_check_model = keep_going

    call check_L_center(s, extras_check_model)
    if (extras_check_model /= keep_going) return

    mass_ratio = M_BH / s% mstar
    mass_ratio_stop = s% x_ctrl(13)
    if (.not. bh_frozen .and. mass_ratio_stop > 0d0 .and. mass_ratio >= mass_ratio_stop) then
        write(*,'(a,f8.5,a,f6.3,a)') '  M_BH/Mstar = ', mass_ratio, ' >= x_ctrl(13) = ', mass_ratio_stop, ': terminating'
        extras_check_model = terminate
        s% termination_code = t_extras_check_model
    end if
end function extras_check_model

integer function extras_finish_step(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_finish_step = keep_going

    if (.not. bh_frozen) then
        call timestep_caps(s, extras_finish_step)
    else if (s% x_ctrl(6) > 0d0) then
        call expand_stop(s, extras_finish_step)
    end if

    if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
end function extras_finish_step

integer function how_many_extra_history_columns(id)
    integer, intent(in) :: id
    how_many_extra_history_columns = 20 + 10
end function how_many_extra_history_columns

subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
    integer, intent(in) :: id, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    integer :: nz, k_mid
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    nz = s% nz
    vals(1:n) = -1d99

    names(1) = 'M_BH';        vals(1) = M_BH / Msun
    names(2) = 'Mdot_in';     vals(2) = Mdot_in / Msun
    names(3) = 'Mdot_BH';     vals(3) = Mdot_BH() / Msun
    names(4) = 'L_Edd';       vals(4) = L_Edd_step / Lsun
    names(5) = 'L_i_over_L_Edd';   vals(5) = s% L_center / L_Edd_step
    names(6) = 'L_conv_base_over_L_Edd'
    names(7) = 'L_sat_over_L_Edd'
    if (L_sat_i_last > 0d0) then
        vals(6) = L_conv_i_last / L_Edd_step
        vals(7) = L_sat_i_last / L_Edd_step
    end if

    names(8)  = 'rho_i';            vals(8)  = s% rho(nz)
    names(9)  = 'T_i';              vals(9)  = s% T(nz)
    names(10) = 'prad_i';           vals(10) = s% prad(nz)
    names(11) = 'pgas_i';           vals(11) = s% pgas(nz)
    names(12) = 'cs_i';             vals(12) = s% csound(nz)
    names(13) = 'gamma1_i';         vals(13) = s% gamma1(nz)
    names(14) = 'nabla_ad_i';       vals(14) = s% grada(nz)
    names(15) = 'L_rad_i';          vals(15) = (s% L(nz) - s% L_conv(nz)) / Lsun
    names(16) = 'L_conv_i';         vals(16) = s% L_conv(nz) / Lsun

    k_mid = minloc(abs(s% r(1:nz) - 0.5d0 * (s% R_center + s% r(1))), 1)
    names(17) = 'pgas_div_prad';    vals(17) = s% pgas(k_mid) / s% prad(k_mid)
    names(18) = 'gamma_ad';         vals(18) = s% gamma1(k_mid)

    names(19) = 'cgrav_factor_center'; vals(19) = s% cgrav(nz) / standard_cgrav
    names(20) = 'cgrav_factor_max';    vals(20) = maxval(s% cgrav(1:nz)) / standard_cgrav

    call other_winds_history_columns(id, 20, names, vals, ierr)
end subroutine data_for_extra_history_columns

integer function how_many_extra_profile_columns(id)
    integer, intent(in) :: id
    how_many_extra_profile_columns = 1 + 6
end function how_many_extra_profile_columns

subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
    integer, intent(in) :: id, n, nz
    character (len=maxlen_profile_column_name) :: names(n)
    real(dp) :: vals(nz,n)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    integer :: k
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    vals(:,:) = 0d0

    names(1) = 'L_sat'
    do k = 1, s% nz
        vals(k,1) = L_sat(s, k) / Lsun
    end do

    call other_winds_profile_columns(id, 1, nz, names, vals, ierr)
end subroutine data_for_extra_profile_columns

end module run_star_extras
