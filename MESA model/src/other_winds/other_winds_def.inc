      ! declarations for other_winds
         
         character (len=strlen) :: wind_scheme
         
         integer, parameter :: x_old_wind = 20

         ! s% xtra(x_time_thermal_eq) holds time the star has been in thermal equilibrium
         ! used to determine the start of ZAMS
         integer, parameter :: x_time_thermal_eq = 21

         ! s% lxtra(lx_pre_ZAMS) true if the star has not already properly settled in the ZAMS
         integer, parameter :: lx_pre_ZAMS = 1

         ! s% xtra(x_old_wind) every step we save the used wind mass loss rate here,
         ! this is used to soften too large changes in wind mass loss rates 
         
         !namelist /other_winds/ &
         !   wind_scheme, &
         !   Kudritzki_scaling_factor, &
         !   Grafener_scaling_factor, &
         !   Stern51_scaling_factor
            
      ! end of declarations for other_winds
