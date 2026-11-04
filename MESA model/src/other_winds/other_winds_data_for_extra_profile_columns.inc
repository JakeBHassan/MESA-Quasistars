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

         names(1) = "t_dyn_test"
         names(2) = "L_excess_test"
         names(3) = "E_excess_test"
         names(4) = "E_excess_final"
         names(5) = "E_diff"
         names(6) = "Mloss"
         do k=1,s% nz
            vals(k,1) = t_dyn_test(k)
            vals(k,2) = L_excess_test(k)
            vals(k,3) = E_excess_test(k)
            vals(k,4) = E_excess_final(k)
            vals(k,5) = E_diff(k)
            vals(k,6) = Mloss(k)

         end do
