subroutine other_winds_history_columns(id, i0, names_out, vals_out, ierr)
    integer, intent(in) :: id, i0
    character (len=maxlen_history_column_name), intent(inout) :: names_out(:)
    real(dp), intent(inout) :: vals_out(:)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    character (len=maxlen_history_column_name) :: names(45)
    real(dp) :: vals(45)
    real(dp) :: t_dyn_test, L_excess_test, E_excess_test
    real(dp) :: mdot_erupt, mdot_dutch
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    include 'other_winds/other_winds_data_for_extra_history_columns.inc'
    names_out(i0+1:i0+10) = names(36:45)
    vals_out(i0+1:i0+10) = vals(36:45)
end subroutine other_winds_history_columns

subroutine other_winds_profile_columns(id, j0, nz, names_out, vals_out, ierr)
    integer, intent(in) :: id, j0, nz
    character (len=maxlen_profile_column_name), intent(inout) :: names_out(:)
    real(dp), intent(inout) :: vals_out(:,:)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    character (len=maxlen_profile_column_name) :: names(6)
    real(dp) :: vals(nz,6)
    integer :: k, i
    integer :: taucritidxarr(1)
    integer :: taucritidx
    real(dp), dimension(nz) :: t_dyn_test, L_excess_test, E_excess_test, masked, masked2, E_diff, E_excess_final, E_binding_final, M_above, Mloss
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    include 'other_winds/other_winds_data_for_extra_profile_columns.inc'
    names_out(j0+1:j0+6) = names(1:6)
    vals_out(1:nz, j0+1:j0+6) = vals(1:nz, 1:6)
end subroutine other_winds_profile_columns
