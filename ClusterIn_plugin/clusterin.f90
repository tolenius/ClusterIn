module clusterin

implicit none

contains

subroutine cluster_dynamics(names_vapor,c_vapor,cs_ref,temp,ipr,t_sim,t_tot,j_total,diameter_approx,&
&    c_inout,&
&    naero,dp_aero_lim,dp_aero,mp_aero,c_aero,pres,&
&    c_evap,nmols_evap,c_coag_molec,c_coag_clust,clust_molec,&
&    c_out_bin,comp_out_bin,c_out_all,clust_out_molec,&
&    c_out,c_bins,&
&    t_ss)

use acdc_simulation_setup, only : solve_ss, sources_and_constants                           ! Logicals for the solution approaches
use acdc_system, only : small_set_mode, n_neutral_monomers, n_charges, nclust_max           ! Simulation mode and numbers of vapors and charging states
use acdc_simulation_setup, only : get_system_size, get_vapor_indices                        ! Parameters related to the simulation system size

use acdc_simulation_setup, only : get_bound_index, get_scav_parameters, get_jbin_parameters, get_ambient_params
use acdc_simulation_setup, only : group_size_bins, nbins

use driver_acdc_J, only : acdc_driver, use_solver                                           ! Driver for the numerical integration

    implicit none
    
    ! Input: ambient conditions and simulation time
    character(len=11), dimension(n_neutral_monomers), intent(in) :: names_vapor             ! Vapor names
    real(kind(1.d0)), intent(inout) :: c_vapor(n_neutral_monomers)                          ! Vapor concentrations (1/m^3)
    real(kind(1.d0)), intent(in) :: cs_ref                                                  ! Reference coagulation sink (when no explicit sink from the aerosol distribution is used) (1/s)
    real(kind(1.d0)), intent(in) :: temp                                                    ! Temperature (K)
    real(kind(1.d0)), intent(in) :: ipr                                                     ! Ion production rate (1/s/m^3)
    real(kind(1.d0)), intent(in) :: t_sim, t_tot                                            ! Simulation time and total accumulated time (s)
    ! Optional input
    real(kind(1.d0)), intent(inout), optional :: c_inout(:)                                 ! Cluster concentration input and output, in case there is need to override the concentrations saved within the routine (1/m^3)
    integer, intent(in), optional :: naero                                                  ! Number of aerosol bins
    real(kind(1.d0)), intent(in), optional :: dp_aero_lim(:), dp_aero(:), mp_aero(:), c_aero(:) ! Aerosol bin diameters (m), masses (kg) and concentrations (1/m^3)
    real(kind(1.d0)), intent(in), optional :: pres                                          ! Ambient pressure (for scavenging kernels) (Pa)
    real(kind(1.d0)), intent(in), optional :: c_evap                                        ! Concentration of particles evaporating back to clusters (1/m^3)
    real(kind(1.d0)), intent(in), optional :: nmols_evap(n_neutral_monomers)                ! Composition of evaporating particles (molec; note: real, not integer)
    
    ! Output: formed particles/time and their acid content
    real(kind(1.d0)), intent(out) :: j_total                                                ! Simulated formation rate (1/s/m^3)
    real(kind(1.d0)), intent(out) :: diameter_approx                                        ! Mass diameter of the formed particles (m)
    ! Optional output
    real(kind(1.d0)), intent(out), optional :: c_coag_molec(:,:)                            ! Concentrations of coagulated molecules on aerosol bins (size: bins, molec) (1/m^3)
    real(kind(1.d0)), intent(out), optional :: c_coag_clust(:,:)                            ! Concentrations of coagulated clusters on aerosol bins (size: bins, clusters) (1/m^3)
    integer, intent(out), optional :: clust_molec(:,:)                                      ! Cluster composition as numbers of vapor molecules (not specifying their possible charge; size: clusters, molec)
    real(kind(1.d0)), intent(out), optional :: c_out_bin(:)                                 ! Concentrations of outgrown clusters to different aerosol bins (size: bins) (1/m^3)
    real(kind(1.d0)), intent(out), optional :: comp_out_bin(:,:)                            ! Composition of outgrown clusters to different aerosol bins (size: bins, molec) (molec; note: real, not integer)
    real(kind(1.d0)), intent(out), optional :: c_out_all(:)                                 ! Concentrations of outgrown individual clusters (size: nclust_out) (1/m^3)
    integer, intent(out), optional :: clust_out_molec(:,:)                                  ! Outgrown cluster composition as numbers of vapor molecules (not specifying their possible charge; size: nclust_out, molec)
    real(kind(1.d0)), intent(out), optional :: c_out                                        ! Concentration of outgrown particles (1/m^3)
    real(kind(1.d0)), intent(out), optional :: c_bins(nbins)                                ! Size bin concentrations (1/m^3)
    real(kind(1.d0)), intent(out), optional :: t_ss                                         ! Time needed to reach steady state (at the given resolution) (s)

    real(kind(1.d0)), save, allocatable :: c(:), c_saved(:)                                 ! Cluster concentrations (1/m^3)
    real(kind(1.d0)), allocatable :: c_coag(:), c_out_clust(:), source(:), c_tmp(:)
    real(kind(1.d0)) :: j_tmp, t_tmp, real_tmp
    real(kind(1.d0)) :: neg_tmp = -1.d-20
    integer, save, allocatable :: fitted(:,:)
    logical, allocatable :: isconst(:)
    integer, save :: neq_syst = 0, nclust_syst = 0, nout_syst(n_charges) = 0                ! Parameters for the simulation system size
    integer, save :: nclust_out = 0, nrange_coag(2) = 0, nrange_out(2) = 0
    real(kind(1.d0)), save :: diameter_max_syst
    integer, save :: n1vapor(n_neutral_monomers) = 0, nmol_vapor(n_neutral_monomers) = 0    ! Indices for vapor monomers
    character(len=11), dimension(n_neutral_monomers), save :: names_vapor_saved
    real(kind(1.d0)) :: diff_nmols_evap(n_neutral_monomers)
    integer, save :: grouped_indices(0:nbins,nclust_max) = 0                                ! Indices of the clusters belonging to each bin
    integer, save :: nclust_per_bin(0:nbins) = 0                                            ! Total number of clusters in each bin

    real(kind(1.d0)) :: j_out(n_charges)                                                    ! Formation rate vector (for neu only, or for neu, neg and pos)
    logical :: int_ok = .true.                                                              ! .false. if integration fails
    real(kind(1.d0)), save :: t_iter = 1.d-16                                               ! Iteration time step (s) for the Euler method
    integer, save :: ipar(4)                                                                ! Parameters for re-calling the monomer settings and rate constants
    logical, save :: firstcall = .true., found_fitted = .false., set_c_coag = .false., set_jbins = .false., set_c_bins = .false.
    integer :: nclust_bound, i, j, k, n

!--------------------------------------------------------------------------------------------------
! Initializations etc. that are needed only at the first call
!--------------------------------------------------------------------------------------------------
    
    if (firstcall) then
        
        ipar = 0
        
        ! Determine the system size and vapor indices
        call get_system_size(neq_syst=neq_syst,nclust_syst=nclust_syst,nout_syst=nout_syst,diameter_max_syst=diameter_max_syst)
        call get_vapor_indices(names_vapor,n1vapor,nmol_vapor)
        names_vapor_saved = names_vapor
        
        ! Initialize the concentrations
        allocate(c(neq_syst))
        allocate(c_saved(neq_syst))
        c = 0.d0
        
        ! Find and save info on possible fitted concentrations
        allocate(source(neq_syst))
        allocate(isconst(neq_syst))
        allocate(fitted(neq_syst,0:neq_syst))
        call sources_and_constants(neq_syst,source,isconst,fitted)
        if (fitted(1,0) .gt. 0) found_fitted = .true.
        
        allocate(c_tmp(neq_syst))
        
        ! Write out brief information on the simulation settings        
        write(*,*)
        
        write(*,*) 'Solving a set of ', neq_syst, ' equations for vapor species ', names_vapor
        
        if (solve_ss) then
            write(*,*) 'Applying the steady-state assumption'
        else
            write(*,*) 'Running simulations for the given time, i.e. no steady-state assumption applied'
        end if
        
        if (.not. solve_ss) then
            write(*,*)
            write(*,*) '--------------------------------------------------------------------------------------------------------'
            
            if (present(c_inout)) then
                write(*,*) 'Input cluster concentrations received'
                write(*,*) 'Warning: You must ensure that the input is CORRECT and UPDATED after each call to cluster dynamics'
            else
                write(*,*) 'No input cluster concentrations received - tracking the concentrations internally'
                write(*,*) 'Warning: This corresponds to modeling the time evolution of ONLY ONE cluster population, and'
                write(*,*) '         you CANNOT call cluster dynamics e.g. in a loop over different grid cells / vertical layers'
                write(*,*) '         - if this is not what you want, use the input cluster concentration option instead'
            end if
        end if
        
        write(*,*) '--------------------------------------------------------------------------------------------------------'
        write(*,*)
        
    end if
    
!--------------------------------------------------------------------------------------------------
! Input to the ACDC driver, updated at every call
!--------------------------------------------------------------------------------------------------
    
    ! Consider the possibility that the vapors come in different order than in previous call
    if (any(names_vapor .ne. names_vapor_saved)) then
        call get_vapor_indices(names_vapor,n1vapor,nmol_vapor)
        names_vapor_saved = names_vapor
    end if
    
    ! Initialize the rate constants etc. at every call because of the varying ambient conditions
    ipar(1) = 0
    ipar(3) = 0
    
    ! Initialize cluster concentrations, if needed
    if (present(c_inout)) then
        c = c_inout
    end if
    
    ! Settings for output for cluster-aerosol coagulation, if needed
    if (present(c_coag_molec) .or. present(c_coag_clust)) then
        
        if (.not. set_c_coag) then
    
            if (solve_ss) then
                write(*,*) 'Cannot output scavenged concentrations when simulating a steady state'
                stop
            end if
            
            call get_system_size(nrange_coag=nrange_coag)
            allocate(c_coag(nclust_syst))
            
            set_c_coag = .true.
            
            ! If cluster-aerosol scavenging transfer is used, also the reduction of vapors due to clustering must be considered
            write(*,*)
            write(*,*) '--------------------------------------------------------------------------------------------------------'
            write(*,*) 'Warning: Outputting scavenged cluster concentrations - ensure that vapors are reduced during clustering:'
            write(*,*) '         - Vapor concentrations must have been set to vary freely in sources_and_constants'
            write(*,*) '         - Cluster equations must have been generated without the --no_eq option'
            write(*,*) '--------------------------------------------------------------------------------------------------------'
            write(*,*)
            
        end if
        
        ! Set the coagulated cluster concentrations to zero at the beginning of the time step to improve accuracy
        c(nrange_coag(1):nrange_coag(2)) = 0.d0
        
    end if
    
    ! Settings for output for size-classified flux out of the cluster regime, if needed
    if (present(c_out_bin) .or. present(c_out_all)) then
        
        if (.not. set_jbins) then
    
            if (solve_ss) then
                write(*,*) 'Not reasonable to output composition-resolved flux out of the cluster regime &
                &    when simulating a steady state'
                stop
            end if
            
            if ((present(c_out_bin) .and. .not. present(comp_out_bin)) .or. &
                 &    (present(c_out_all) .and. .not. present(clust_out_molec))) then
                write(*,*) 'Not reasonable to output composition-resolved flux out of the cluster regime &
                &    without the respective compositions'
                stop
            end if
            
            call get_system_size(nclust_out=nclust_out,nrange_out=nrange_out)
            allocate(c_out_clust(nclust_out))
            
            set_jbins = .true.
            
        end if
        
        ! Set the outgrown cluster concentrations to zero at the beginning of the time step to improve accuracy
        c(nrange_out(1):nrange_out(2)) = 0.d0
        
    end if
    
    ! Set the outgrown cluster concentrations to zero at the beginning of the time step to improve accuracy
    c(nout_syst) = 0.d0
    
    ! Save the initial situation
    c_saved = c
    
    ! Set the vapor concentrations
    c(n1vapor) = c_vapor
    
    if (present(c_out)) then
        c_out = 0.d0
    end if

!--------------------------------------------------------------------------------------------------
! Optional: Input for aerosol distribution and other parameters for including cluster-aerosol interactions
!--------------------------------------------------------------------------------------------------
    
    ! Update the aerosol parameters in the routine that stores them
    if (present(naero)) call get_ambient_params(naero_in=naero)
    if (present(dp_aero_lim)) call get_ambient_params(dp_aero_lim_in=dp_aero_lim)
    if (present(dp_aero)) call get_ambient_params(dp_aero_in=dp_aero)
    if (present(mp_aero)) call get_ambient_params(mp_aero_in=mp_aero)
    if (present(c_aero)) call get_ambient_params(c_aero_in=c_aero)
    if (present(pres)) call get_ambient_params(pres_in=pres)
    
!--------------------------------------------------------------------------------------------------
! Optional: Input for evaporation of larger aerosols to the cluster regime
!--------------------------------------------------------------------------------------------------
    
    ! Insert evaporated larger particles, if needed
    if (present(c_evap)) then
        if (c_evap .gt. 0.d0) then
            
            if (solve_ss) then
                write(*,*) 'Cannot input evaporated aerosols when simulating a steady state'
                stop
            end if
            
            if (present(nmols_evap)) then
                
                ! Insert the evaporated particles to the cluster regime
                call get_bound_index(nmol_vapor,nmols_evap,nclust_bound,diff_nmols_evap)
                c(nclust_bound) = c(nclust_bound) + c_evap
                ! Account for possible differences in the exact composition by adjusting the vapor concentrations 
                c_vapor = c_vapor + c_evap*diff_nmols_evap
                
                ! Subtract c_evap from the outgrown concentration for consistency in standalone tests
                ! c(nout_syst(1)) = c(nout_syst(1)) - c_evap ! Assume neutral channel
                ! if (present(c_out)) then
                    ! c_out = c_out - c_evap
                ! end if
                
            else
                write(*,*) 'If c_evap is used, also nmols_evap must be given to cluster_dynamics'
                stop
            end if
            
        end if
    end if
    
!--------------------------------------------------------------------------------------------------
! Call the ACDC driver and get the standard output
!--------------------------------------------------------------------------------------------------
    
    ! Run the simulation to the steady state / for the given time
    t_tmp = t_sim
    call acdc_driver(neq_syst,nclust_syst,nout_syst,c,cs_ref,temp,ipr,t_tmp,t_tot,t_iter,ipar,int_ok,j_out)
    
    ! Get the formation rate
    j_total = sum(j_out)
    !write(*,*) 'Neutral, neg. and pos. J: ',j_out*1.d-6,' cm^-3 s^-1'
    
    ! Update the vapor concentrations as clusters may be assumed to act as a sink for vapors
    ! Note that if some concentrations are fitted (e.g. [H2SO4] -> sum([H2SO4*base])), the returned value is
    ! the sum instead of the actual monomer value
    if (.not. found_fitted) then
        c_vapor = c(n1vapor)
    else
        ! If some vapors have fitted concentrations, store the sum values in a temporary concentration array
        c_tmp = c
        do n = 1,fitted(1,0)
            k = fitted(n+1,0)
            c_tmp(k) = sum(c((/k, fitted(n+1,2:fitted(n+1,1)+1)/)))
        end do
        c_vapor = c_tmp(n1vapor)
    end if
    
    ! The first outgrown size is ~equal to the largest simulated size
    diameter_approx = diameter_max_syst
    
!--------------------------------------------------------------------------------------------------
! Optional: Output for cluster-aerosol coagulation
!--------------------------------------------------------------------------------------------------
    
    ! Determine the numbers of clusters or molecules transferred to each aerosol bin through cluster scavenging
    if (present(c_coag_molec) .or. present(c_coag_clust)) then
        
        ! Clusters scavenged during the simulation time step
        c_coag = c(nrange_coag(1):nrange_coag(2)) - c_saved(nrange_coag(1):nrange_coag(2))
        
        real_tmp = minval(c_coag)
        if (real_tmp .lt. 0.d0) then
            if (real_tmp .lt. neg_tmp) then
                write(*,*) 'Warning: negative coagulated concentrations, lowest value = ',real_tmp*1.d-6,' cm^-3, something wrong?'
            end if
            where(c_coag .lt. 0.d0) c_coag = 0.d0
        end if
        
        if (present(c_coag_molec)) then
            ! Simplified approach (usually the default): Return the concentration of each vapor species
            call get_scav_parameters(names_vapor=names_vapor,temperature=temp,c_coag=c_coag,&
            &    molec_per_aero=c_coag_molec)
        else
            ! Accurate approach (for model testing): Return the concentration of each cluster plus the cluster compositions
            call get_scav_parameters(names_vapor=names_vapor,temperature=temp,c_coag=c_coag,&
            &    clust_per_aero=c_coag_clust,nmols_nocharge_clust=clust_molec)
            
        end if
        
    end if

!--------------------------------------------------------------------------------------------------
! Optional: Output for size-classified flux out of the cluster regime
!--------------------------------------------------------------------------------------------------
    
    if (present(c_out_bin) .or. present(c_out_all)) then
        
        ! Clusters grown out to different bins during the simulation time step
        c_out_clust = c(nrange_out(1):nrange_out(2)) - c_saved(nrange_out(1):nrange_out(2))
        
        real_tmp = minval(c_out_clust)
        if (real_tmp .lt. 0.d0) then
            if (real_tmp .lt. neg_tmp) then
                write(*,*) 'Warning: negative outgrown concentrations, lowest value = ',real_tmp*1.d-6,' cm^-3, something wrong?'
            end if
            where(c_out_clust .lt. 0.d0) c_out_clust = 0.d0
        end if
        
        if (present(c_out_all)) c_out_all = c_out_clust
        
        if (present(clust_out_molec)) then
            call get_jbin_parameters(c_out_clust,c_out_bin,comp_out_bin,nmols_nocharge_out=clust_out_molec,names_vapor=names_vapor)
        else
            call get_jbin_parameters(c_out_clust,c_out_bin,comp_out_bin,names_vapor=names_vapor)
        end if
        
        j_tmp=sum(c_out_bin)/t_sim
        
        if ((j_total .gt. 1.d-6) .or. (j_tmp .gt. 1.d-6)) then
            if (abs((j_total-j_tmp)/max(j_total,j_tmp)) .gt. 1.d-2) then
                write(*,*) 'Warning: the sum of binned flux out does not equal the overall flux, something wrong?'
                write(*,*) j_tmp*1.d-6,' cm^-3 s^-1 vs ',j_total*1.d-6,' cm^-3 s^-1'
            end if
        end if
        
        ! For consistency, set the outputted J to the sum of the binned J (as there may be very minor differences due to the numerical integration)
        j_total=j_tmp
        
    end if
    
    ! See that the formation rate is ok
    if (j_total .lt. 0.d0) then
        if (j_total .lt. neg_tmp) then
            write(*,*) 'Warning: negative J = ',j_total*1.d-6,' cm^-3 s^-1, something wrong?'
        end if
        j_total = 0.d0
    end if
    
!--------------------------------------------------------------------------------------------------
! Optional: Output for cluster concentrations
!--------------------------------------------------------------------------------------------------
    
    ! Simulated clusters
    if (present(c_inout)) then
        c_inout = c
    end if
    
    ! Outgrown clusters
    if (present(c_out)) then
        ! Cumulative c_out
        !c_out = sum(c(nout_syst)) ! Arbitrary for solve_ss = .true.
        ! Clusters grown out during the simulation time step
        c_out = c_out + j_total*t_sim
    end if
    
    ! Simulated clusters grouped into size bins
    if (present(c_bins)) then
        if (.not. small_set_mode) then
            
            if (.not. set_c_bins) then
                call group_size_bins(grouped_indices,nclust_per_bin)
                set_c_bins = .true.
            end if
            
            ! Concentration in bins
            c_bins = 0.d0
            do i = 1,nbins
                if (nclust_per_bin(i) .ne. 0) then
                    c_bins(i) = sum(c(grouped_indices(i,1:nclust_per_bin(i))))
                end if
            end do
            
        end if
    end if
    
!--------------------------------------------------------------------------------------------------
! Optional: Output for steady-state relaxation time scale
!--------------------------------------------------------------------------------------------------
    
    if (present(t_ss)) then
        t_ss = t_tmp
    end if

!--------------------------------------------------------------------------------------------------
! Info after the processing that is needed only at the first call
!--------------------------------------------------------------------------------------------------
    
    if (firstcall) then
    
        firstcall = .false.
        
        write(*,*)
        write(*,*) '--------------------------------------------------------------------------------------------------------'
        if (use_solver) then
            write(*,*) 'VODE solver chosen for integrating the cluster equations'
        else
            write(*,*) 'Large cluster system: a simpler Euler method chosen instead of VODE solver'
        end if
        write(*,*) '--------------------------------------------------------------------------------------------------------'
        write(*,*)
        
    end if
    
end subroutine cluster_dynamics

end module clusterin
