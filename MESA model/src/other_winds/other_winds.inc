            pure function integrate(x, y) result(r)
              real(dp), intent(in)  :: x(:)
              real(dp), intent(in)  :: y(size(x))
              real(dp)              :: r

              associate(n => size(x))
                r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
              end associate
            end function

            subroutine get_mloss_forwind(id,ierr, testsize,mdot_erupt)
                 type (star_info), pointer :: s
                 integer, intent(in) :: testsize, id
                 integer, intent(out) :: ierr
                 integer :: taucritidxarr(1)
                 integer :: taucritidx
                 integer :: i,  mlost_loc
                 real(dp) :: mlost_final, tdyn_final, masslossres, tdyn_final_alt, masslossres_alt, mdot_erupt
                 real(dp), dimension(testsize) ::  t_dyn_test, L_excess_test, E_excess_test, E_excess_pos, masked, masked2, E_diff, E_excess_final, E_binding_final, E_diff_masked, M_above, mdot_loss, Mloss
                 include 'formats'

                 ierr = 0
                 call star_ptr(id, s, ierr)
                 if (ierr /= 0) return

                 taucritidxarr = minloc(abs(s% tau(1:s% nz)-(clight/(s% csound(1:s% nz)))))
                 taucritidx = taucritidxarr(1)

                 do i  =1,s% nz
                    t_dyn_test(i) = 0d0
                    L_excess_test(i) = 0d0
                    E_excess_test(i) = 0d0
                 end do

                 t_dyn_test = ((s% r(1:s% nz))**3/(6.674e-8*s% m(1:s% nz)))**(1.0/2.0)
                 L_excess_test = (s% L(1:s% nz)) - s% L_conv(1:s% nz) - (pi4*clight*s% cgrav(1:s% nz)*s% m_grav(1:s% nz)/(s% opacity(1:s% nz)))
                 E_excess_test = L_excess_test*t_dyn_test

                 where (E_excess_test < 0)
                     masked2 = 0
                elsewhere
                     masked2 = E_excess_test
                end where

                where ((s% tau(1:s% nz)) > (clight/(s% csound(1:s% nz))))
                     masked = 0
                elsewhere
                     masked = masked2
                end where

                do i  =1,s% nz-1
                     E_excess_final(s% nz+1-i) = -1*integrate(s% m(s% nz-i:s% nz), masked(s% nz-i:s% nz))/ABS(s% m(s% nz+1-i)-s% m(taucritidx))
                end do

                do i = 1, s% nz-1
                     E_binding_final(i) = SUM(6.674e-8*s% m(1:i)*s% dm(1:i)/s% r(1:i))
                     M_above(i) = (s% m(1) - s% m(i))/1.9884098706980504d33
                end do

                E_excess_final(1) = E_excess_final(2)
                E_binding_final(s% nz) = E_binding_final(s% nz-1)
                M_above(s% nz) = (s% m(1) - s% m(s% nz))/1.9884098706980504d33

                E_diff = E_excess_final-E_binding_final

                where (E_diff > 0)
                    Mloss = M_above
                elsewhere
                    Mloss = 0d0
                end where

                Mloss(taucritidx) = 0d0

                mlost_loc = MAXLOC(Mloss, DIM=1)
                mlost_final = MAXVAL(Mloss)

                tdyn_final = t_dyn_test(mlost_loc)*3.17098d-8

                masslossres = mlost_final/tdyn_final

                mdot_erupt = masslossres

                s% xtra(26) = M_above(taucritidx)
                s% xtra(27) = mlost_final
                s% xtra(28) = tdyn_final
                s% xtra(29) = mlost_loc

            end subroutine get_mloss_forwind

            subroutine erupt_other_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, wind, ierr)
               use star_def
               integer, intent(in) :: id
               real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z
               real(dp), intent(out) :: wind
               integer, intent(out) :: ierr

               integer :: h1, he4
               real(dp) :: w, Xs, Ys, Z_div_Z_solar, Teff_jump, alfa, L1, M1, R1, T1, &
                  vink_wind, nieu_wind, hamann_wind, highT_w, lowT_w, Twindow, T_high, T_low, surface_h1, w1, w2, T_decin
               real(dp), parameter :: Zsolar = 0.019d0

               type (star_info), pointer :: s
               ierr = 0
               call star_ptr(id, s, ierr)
               if (ierr /= 0) return

               L1 = Lsurf
               M1 = Msurf
               R1 = Rsurf
               T1 = Tsurf

               h1 = s% net_iso(ih1)
               he4 = s% net_iso(ihe4)
               Xs = s% xa(h1,1)
               Ys = s% xa(he4,1)
               surface_h1 = s% xa(h1,1)

               T_high = 11000
               T_low = 10000
               T_decin = 4000
               if (s% Dutch_scaling_factor == 0) then
                  wind = 0
               else if (T1 <= T_low) then
                   if (T1 > T_decin) then
                     call eval_lowT_Dutch(wind)
                   else
                     call eval_decin_wind(wind)
                 end if
               else if (T1 >= T_high) then
                  call eval_highT_Dutch(wind)
               else
                  call eval_lowT_Dutch(w1)
                  call eval_highT_Dutch(w2)
                  alfa = (T1 - T_low)/(T_high - T_low)
                  wind = (1-alfa)*w1 + alfa*w2
               end if

               s% xtra(23) = s% Dutch_scaling_factor*wind

               w = 0

               call get_mloss_forwind(id, ierr, s% nz, w)

               s% xtra(25) = w * s% x_ctrl(15)

               print*, 'Mdot erupt', s% xtra(25)

               w = w * s% x_ctrl(15)

               wind = w + s% Dutch_scaling_factor*wind

               s% xtra(24) = wind

               contains

               subroutine eval_highT_Dutch(w)
                  real(dp), intent(out) :: w
                  include 'formats'
                  if (surface_h1 < 0.4d0) then
                     w = 1d-11 * pow(L1/Lsun,1.29d0) * pow(Y,1.7d0) * sqrt(Z)
                  else
                     call eval_Vink_wind(w)
                  end if
               end subroutine eval_highT_Dutch

               subroutine eval_lowT_Dutch(w)
                  real(dp), intent(out) :: w
                  include 'formats'
                  if (s% Dutch_wind_lowT_scheme == 'de Jager') then
                     call eval_de_Jager_wind(w)
                  else
                     write(*,*) 'unknown value for Dutch_wind_lowT_scheme ' // &
                        trim(s% Dutch_wind_lowT_scheme)
                     w = 0
                  end if
               end subroutine eval_lowT_Dutch

               subroutine eval_Vink_wind(w)
                  real(dp), intent(inout) :: w
                  real(dp) :: alfa, w1, w2, Teff_jump, logMdot, dT, vinf_div_vesc

                  if (T1 > 27500d0) then
                     alfa = 1
                  else if (T1 < 22500d0) then
                     alfa = 0
                  else
                     Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z/Zsolar)))
                     dT = 100d0
                     if (T1 > Teff_jump + dT) then
                        alfa = 1
                     else if (T1 < Teff_jump - dT) then
                        alfa = 0
                     else
                        alfa = 0.5d0*(T1 - (Teff_jump - dT)) / dT
                     end if
                  end if

                  if (alfa > 0) then
                     vinf_div_vesc = 2.6d0
                     vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0)
                     logMdot = &
                        - 6.697d0 &
                        + 2.194d0*log10(L1/Lsun/1d5) &
                        - 1.313d0*log10(M1/Msun/30d0) &
                        - 1.226d0*log10(vinf_div_vesc/2d0) &
                        + 0.933d0*log10(T1/4d4) &
                        - 10.92d0*pow2(log10(T1/4d4)) &
                        + 0.85d0*log10(Z/Zsolar)
                     w1 = exp10(logMdot)
                  else
                     w1 = 0
                  end if

                  if (alfa < 1) then
                     vinf_div_vesc = 1.3d0
                     vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0)
                     logMdot = &
                        - 6.688d0 &
                        + 2.210d0*log10(L1/Lsun/1d5) &
                        - 1.339d0*log10(M1/Msun/30d0) &
                        - 1.601d0*log10(vinf_div_vesc/2d0) &
                        + 1.07d0*log10(T1/2d4) &
                        + 0.85d0*log10(Z/Zsolar)
                     w2 = exp10(logMdot)
                  else
                     w2 = 0
                  end if

                  w = alfa*w1 + (1 - alfa)*w2

               end subroutine eval_Vink_wind

               subroutine eval_de_Jager_wind(w)
                  real(dp), intent(out) :: w
                  real(dp) :: log10w
                  include 'formats'
                  log10w = 1.769d0*log10(L1/Lsun) - 1.676d0*log10(T1) - 8.158d0
                  w = exp10(log10w)
               end subroutine eval_de_Jager_wind

               subroutine eval_decin_wind(w)
                  real(dp), intent(out) :: w
                  real(dp) :: log10w, term1, term2, term3
                  include 'formats'
                  term1 = 1.77d0+0.49d0
                  term2 = -1.68d0 * (s% initial_mass/10d0)
                  term3 = 3.5d0*log10(L1/(1d5*Lsun))

                  log10w = term1+term2+term3

                  w = exp10(log10w)*(1d-5)
               end subroutine eval_decin_wind

            end subroutine erupt_other_wind

      subroutine brott_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z
         real(dp), intent(out) :: w
         integer, intent(out) :: ierr

         integer :: h1, he4
         real(dp) :: Xs, Ys, Z_div_Z_solar, Teff_jump, alfa, L1, M1, R1, T1, &
            vink_wind, nieu_wind, hamann_wind, highT_w, lowT_w, Twindow
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         L1 = Lsurf
         M1 = Msurf
         R1 = Rsurf
         T1 = Tsurf

         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         Xs = s% xa(h1,1)
         Ys = s% xa(he4,1)
         Z_div_Z_solar = s% kap_rq% Zbase/0.0142d0
         Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z_div_Z_solar)))

         vink_wind = 0d0
         nieu_wind = 0d0
         hamann_wind = 0d0
         w = 0

         call eval_Vink_wind(vink_wind)
         call eval_Nieuwenhuijzen_wind(nieu_wind)
         call eval_Hamann_wind(hamann_wind)

         hamann_wind = hamann_wind/10d0

         lowT_w = max(vink_wind, nieu_wind)

         alfa = 0d0
         if (Xs > 0.7d0) then
            alfa = 1d0
         else if (Xs > 0.4d0 .and. Xs < 0.7d0) then
            alfa = (Xs - 0.4d0)/0.3d0
         end if
         highT_w = alfa * vink_wind + (1d0-alfa) * hamann_wind

         Twindow = Teff_jump*0.10d0
         alfa = 0d0
         if (T1 < Teff_jump - Twindow/2d0) then
            alfa = 1d0
         else if (T1 > Teff_jump - Twindow/2d0 .and. T1 < Teff_jump + Twindow/2d0) then
            alfa = ((Teff_jump + Twindow/2d0)-T1)/Twindow
         end if
         w = alfa * lowT_w + (1d0-alfa) * highT_w

         if (s% xtra(x_old_wind) /= 0) then
            if(abs(w) > abs(s% xtra(x_old_wind))*1.05) then
               w = s% xtra(x_old_wind)*1.05
            else if(abs(w) < abs(s% xtra(x_old_wind))*0.95) then
               w = s% xtra(x_old_wind)*0.95
            end if
         end if
         s% xtra(x_old_wind) = w

         ierr = 0

         contains

         subroutine eval_Vink_wind(w)
            real(dp), intent(inout) :: w
            real(dp) :: alfa, w1, w2, logMdot, dT, vinf_div_vesc

            if (T1 > 27500d0) then
               alfa = 1
            else if (T1 < 22500d0) then
               alfa = 0
            else
               dT = 100d0
               if (T1 > Teff_jump + dT) then
                  alfa = 1
               else if (T1 < Teff_jump - dT) then
                  alfa = 0
               else
                  alfa = (T1 - (Teff_jump - dT)) / (2*dT)
               end if
            end if

            if (alfa > 0) then
               vinf_div_vesc = 2.6d0
               vinf_div_vesc = vinf_div_vesc*pow(Z_div_Z_solar,0.13d0)
               logMdot = &
                  - 6.697d0 &
                  + 2.194d0*log10(L1/Lsun/1d5) &
                  - 1.313d0*log10(M1/Msun/30) &
                  - 1.226d0*log10(vinf_div_vesc/2d0) &
                  + 0.933d0*log10(T1/4d4) &
                  - 10.92d0*pow2(log10(T1/4d4)) &
                  + 0.85d0*log10(Z_div_Z_solar)
               w1 = 10**(logMdot)
            else
               w1 = 0
            end if

            if (alfa < 1) then
               vinf_div_vesc = 1.3d0
               vinf_div_vesc = vinf_div_vesc*pow(Z_div_Z_solar,0.13d0)
               logMdot = &
                  - 6.688d0 &
                  + 2.210d0*log10(L1/Lsun/1d5) &
                  - 1.339d0*log10(M1/Msun/30) &
                  - 1.601d0*log10(vinf_div_vesc/2d0) &
                  + 1.07d0*log10(T1/2d4) &
                  + 0.85d0*log10(Z_div_Z_solar)
               w2 = 10**(logMdot)
            else
               w2 = 0
            end if

            w = alfa*w1 + (1 - alfa)*w2

         end subroutine eval_Vink_wind

         subroutine eval_Nieuwenhuijzen_wind(w)
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -14.02d0 &
                     +1.24d0*log10(L1/Lsun) &
                     +0.16d0*log10(M1/Msun) &
                     +0.81d0*log10(R1/Rsun)
            w = 10**(log10w)
         end subroutine eval_Nieuwenhuijzen_wind

         subroutine eval_Hamann_wind(w)
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -11.95d0 &
                     +1.5d0*log10(L1/Lsun) &
                     -2.85d0*Xs &
                     + 0.85d0*log10(Z_div_Z_solar)
            w = 10**(log10w)
         end subroutine eval_Hamann_wind

      end subroutine brott_wind
