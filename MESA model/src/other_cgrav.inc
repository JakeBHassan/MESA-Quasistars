subroutine my_other_cgrav(id, ierr)
    use utils_lib, only: is_bad
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type(star_info), pointer :: s
    integer :: k
    integer, parameter :: tov_ramp_start = 100, tov_ramp_models = 50
    real(dp) :: r, G, m, rho, P, u, f1, f2, f3, ramp, taper
    include 'formats'
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    ramp = min(1d0, max(0d0, dble(s% model_number - s% init_model_number - tov_ramp_start) / tov_ramp_models))
    if (s% doing_relax) ramp = 0d0
    do k=1,s% nz
        r = s% r(k)
        G = standard_cgrav
        m = s% m(k)
        if (k > 1) then
            rho = 0.5d0*(s% rho(k-1) + s% rho(k))
            P = 0.5d0*(s% Peos(k-1) + s% Peos(k))
            u = 0.5d0*(s% energy(k-1) + s% energy(k))
            f1 = 1d0 + (u + P/rho)/clight**2
            f2 = 1d0 + pi4*pow3(r)*P/(m*clight**2)
        else
            f1 = 1d0
            f2 = 1d0
        end if
        f3 = 1d0 - 2d0*G*m/(r*clight**2)

        taper = min(1d0, (2d0*G*m/(s% r_start(k)*clight**2)/1d-4)**2)
        s% cgrav(k) = G*(1d0 + ramp*taper*(f1*f2/f3 - 1d0))
        if (is_bad(s% cgrav(k))) then
            write(*,2) 's% cgrav(k)', k, s% cgrav(k)
            write(*,2) 'f1', k, f1
            write(*,2) 'f2', k, f2
            write(*,2) 'f3', k, f3
            write(*,2) 'r', k, r
            write(*,2) 'G', k, G
            write(*,2) 'm', k, m
            if (k > 1) then
                write(*,2) 'P', k, P
                write(*,2) 'rho', k, rho
            end if
            stop 'my_other_cgrav'
        end if
    end do
end subroutine my_other_cgrav
