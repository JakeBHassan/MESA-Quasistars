        t_dyn_test = 0d0
        L_excess_test = 0d0
        E_excess_test = 0d0
        
        mdot_erupt = 0d0
        mdot_dutch = 0d0
        
        call get_mloss_forwind(id, ierr, s% nz, mdot_erupt)
        
        t_dyn_test = ((s% r(1))**3/(6.674e-8*s% m(1)))**(1.0/2.0)
        L_excess_test = (s% L(1)) - (pi4*clight*s% cgrav(1)*s% m_grav(1)/(s% opacity(1)))
        E_excess_test = L_excess_test*t_dyn_test
        
        names(36) = 't_dyn_test'
        vals(36) = t_dyn_test
        
        names(37) = 'L_excess_test'
        vals(37) = L_excess_test
        
        names(38) = 'E_excess_test'
        vals(38) = E_excess_test
        
        names(39) = 'mdot_erupt'
        vals(39) = s% xtra(25)
        
        names(40) = 'mdot_dutch'
        vals(40) = s% xtra(23)
        
        names(41) = 'mdot_test_total'
        vals(41) = s% xtra(24)
        
        names(42) = 'Mabove_crit'
        vals(42) = s% xtra(26)
        
        names(43) = 'mlost_final'
        vals(43) = s% xtra(27)
        
        names(44) = 'tdyn_final'
        vals(44) = s% xtra(28)
        
        names(45) = 'mlost_loc'
        vals(45) = s% xtra(29)