        !Find the index at which critical optical depth occurs
        taucritidxarr = minloc(abs(s% tau(1:s% nz)-(clight/(s% csound(1:s% nz)))))
        taucritidx = taucritidxarr(1)
         
        !Initialize arrays
        do i  =1,s% nz
        t_dyn_test(i) = 0d0
        L_excess_test(i) = 0d0
        E_excess_test(i) = 0d0
        end do
         
        !Calculate dynamical time, excess Luminosity, and excess Energy. 
        t_dyn_test = ((s% r(1:s% nz))**3/(6.674e-8*s% m(1:s% nz)))**(1.0/2.0)
        L_excess_test = (s% L(1:s% nz)) - s% L_conv(1:s% nz) - (pi4*clight*s% cgrav(1:s% nz)*s% m_grav(1:s% nz)/(s% opacity(1:s% nz)))
        E_excess_test = L_excess_test*t_dyn_test
         
        !Mask out negative excess Energy
        where (E_excess_test < 0)
            masked2 = 0
        elsewhere
            masked2 = E_excess_test
        end where
        
        !Mask out where tau < critical tau
        where ((s% tau(1:s% nz)) > (clight/(s% csound(1:s% nz))))
            masked = 0 
        elsewhere 
            masked = masked2
        end where
        
        !Calculate total excess energy from surface to critical tau 
        do i  =1,s% nz-1
            E_excess_final(s% nz+1-i) = -1*integrate(s% m(s% nz-i:s% nz), masked(s% nz-i:s% nz))/ABS(s% m(s% nz+1-i)-s% m(taucritidx)) !mass method
        end do
        
        !Find binding energy, mass above
        do i = 1, s% nz-1
            E_binding_final(i) = SUM(6.674e-8*s% m(1:i)*s% dm(1:i)/s% r(1:i)) !still right way, surface first
            M_above(i) = (s% m(1) - s% m(i))/1.9884098706980504d33
        end do
        
        !Patch up first element of arrays
        E_excess_final(1) = E_excess_final(2)
        E_binding_final(s% nz) = E_binding_final(s% nz-1)
        M_above(s% nz) = (s% m(1) - s% m(s% nz))/1.9884098706980504d33
        
        !Define energy difference between excess (due to local supereddington L) and binding
        E_diff = E_excess_final-E_binding_final

        !Total lost mass: only where E_diff > 0
        where (E_diff > 0)
            Mloss = M_above
        elsewhere
            Mloss = 0d0
        end where
        
        Mloss(taucritidx) = 0d0
         
         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.
         !names(1) = "Lrad_div_Leddd"
         !names(2) = "Lrad_div_Leddd2"
         !names(3) = "Ledd_opacity"
         !names(4) = "Ledd_electron"
         !names(5) = "Gamma_opacity"
         !names(6) = "Gamma_electron"
         !do k=1,s% nz
         !vals(k,1) = 4d0*crad/3d0*pow4(s% T(k))/s% Peos(k)*s%gradT(k)
         !vals(k,2) = 4d0*crad/3d0*pow4(s% T(k))/s% Peos(k)*s%mlt_gradT(k)
         !vals(k,3) = pi4*clight*s% cgrav(k)*s% m_grav(k)/(s% opacity(k)*Lsun)
         !vals(k,4) = pi4*clight*s% cgrav(k)*s% m_grav(k)/(0.2d0*(1+s% x(k))*Lsun)
         !vals(k,5) = s% L(k) / vals(k,3)
         !vals(k,6) = s% L(k) / vals(k,4)
         !! vals(k,7) = pi4*clight*s% cgrav(k)*s% m_grav(k)/(*Lsun)
         
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
            ! vals(k,7) = pi4*clight*s% cgrav(k)*s% m_grav(k)/(*Lsun)

         end do