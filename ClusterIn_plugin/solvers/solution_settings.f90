module solution_settings

implicit none

    ! Maximum number of equations for using a solver; otherwise a simpler Euler method is used
    integer, parameter :: max_eq_solver = 5000             ! (Set to a high number just to favor the solver; should probably be decreased)
    
    ! Relative tolerance in the solver or when assessing the Euler solutions (see driver)
    real(kind(1.d0)), parameter :: rtol_solver = 1.d-5
    !real(kind(1.d0)), parameter :: rtol_solver = 1.d-3    ! 1.d-3 is the default tolerance e.g. in Matlab
    
    ! Settings for a solver
    real(kind(1.d0)), parameter :: atol_solver = 1.d-20    ! Absolute tolerance in the solver
    
    ! Settings for Euler integration
    real(kind(1.d0)), parameter :: chmax = 1.d-2           ! Threshold relative change in the concentration for
                                                           ! determining if the time step should be decreased
    real(kind(1.d0)), parameter :: chtol = 1.d-20          ! Lowest conc. (m^-3) to consider when comparing the changes
    
    ! Lowest negative concentration accepted
    real(kind(1.d0)), parameter :: negtol = -1.d-20        ! (m^-3)
    
    ! Criteria for monitoring the convergence to a steady state
    real(kind(1.d0)), parameter :: sstol = 1.d-5           ! Maximum relative change in the concentrations
    real(kind(1.d0)), parameter :: sstime_res = 6.d2       ! Minimum time interval for checking the changes (s)
    !real(kind(1.d0)), parameter :: sstime_res = 60.d0
    real(kind(1.d0)), parameter :: sstime_tot = 1.2d3      ! Minimum total simulation time required for a steady state (s)
    ! Maximum time allowed for trying to reach steady state (set in order to not get stuck at extremely slowly converging conditions)
    real(kind(1.d0)), parameter :: sstime_max = 1.d8       ! Maximum total simulation time after which to give up (s)
    logical, parameter :: l_sstime = .false.               ! Monitor also the time needed to reach steady state (at the resolution of sstime_res) (s)


end module solution_settings
