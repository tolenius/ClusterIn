program run_clusterin_example

use acdc_aerosol_parameters, only : ctot_aero, dmed_aero, sigma_aero    ! Aerosol parameters for testing the dynamic coupling
use acdc_aerosol_parameters, only : dp_lim_aero, lognormal_distr        ! (will in practice be given by the aerosol dynamics model)
use acdc_aerosol_parameters, only : naero => naero_fixed

use clusterin, only : cluster_dynamics                                  ! The main module for cluster input / output
use acdc_simulation_setup, only : get_system_size                       ! For testing the dynamic features

    implicit none
    
    ! All variables in SI units
    
    integer :: n_clustering_vapors                                  ! Number of clustering vapors
    character(len=11), dimension(:), allocatable :: names_vapor     ! Vapor names for clustering species
    real(kind(1.d0)), allocatable :: c_vapor(:)                     ! Vapor concentrations for clustering species (m^-3)
    real(kind(1.d0)) :: temp                                        ! Temperature (K)
    real(kind(1.d0)) :: ipr                                         ! Ion production rate (m^-3 s^-1)
    real(kind(1.d0)) :: c_out                                       ! Total concentration of formed particles during the time step (m^-3)
    
    ! Parameters for the "traditional" set-up (dummies if the "dynamic" set-up is used)
    
    real(kind(1.d0)) :: cs_ref                                      ! Coagulation sink (when not determined explicitly from the aerosol distribution) (s^-1)
    real(kind(1.d0)) :: j_acdc                                      ! Simulated total formation rate (m^-3 s^-1)
    real(kind(1.d0)) :: diameter_acdc                               ! Approximate mass diameter of the formed particles (m)
    
    ! Parameters for the "dynamic" set-up
    
    ! Mass transfer from clusters to aerosols through cluster scavenging
    character(len=15) :: coag_test_name = ''                        ! Name of the scavenging mass transfer test (either molecules or clusters transferred)
    integer :: nclust_syst                                          ! Number of clusters, possibly needed for testing the dynamic features
    real(kind(1.d0)) :: dp_aero_lim(naero+1), dp_aero(naero), c_aero(naero)  ! Aerosol parameters, needed for the dynamic features
    real(kind(1.d0)), allocatable  :: c_coag_aero(:,:)              ! Concentrations of coagulated molecules or clusters on aerosol bins (1/m^3)
    integer, allocatable  :: clust_molec(:,:)                       ! Cluster composition as numbers of vapor molecules (not specifying their possible charge)
    
    ! Formed particles that may be of different sizes (due to formation through different cluster-molecule and cluster-cluster collisions)
    real(kind(1.d0)), allocatable  :: c_out_bin(:)                  ! Concentrations of outgrown clusters to different aerosol bins (1/m^3)
    real(kind(1.d0)), allocatable :: comp_out_bin(:,:)              ! Composition of outgrown clusters to different aerosol bins (molec; note: real, not integer)
    
    ! Evaporation of larger particles back to the cluster regime
    real(kind(1.d0)) :: c_evap                                      ! Concentration of particles evaporating back to cluster regime (m^-3)
    real(kind(1.d0)), allocatable :: nmols_evap(:)                  ! Composition of evaporating particles (molec; note: real, not integer)
    
    ! Misc.
    
    real(kind(1.d0)) :: time_res, t                                 ! Simulation times (s)
    
    ! Some helper variables
    logical :: l_dynamic
    character(len=11), dimension(:), allocatable :: str_vapor
    character(len=20) :: fmt_es = '(*(es15.2e4,1x))'
    character(len=200) :: char_tmp
    integer :: i, j, k, n, it
    
    
    ! Input
    
    ! Simulation conditions
    
    ! Let's say that the cluster system includes 2 vapor species
    n_clustering_vapors = 2
    allocate(names_vapor(n_clustering_vapors))
    allocate(c_vapor(n_clustering_vapors))
    names_vapor(1)(:) = 'A'
	names_vapor(2)(:) = 'N'
	c_vapor = (/2.d6, 1.d9/)*1.d6
    
    temp = 280.d0
    ipr = 3.d0*1.d6 ! When ions are not included, IPR is only a dummy
    cs_ref = 1.d-3  ! For the dynamic set-up, CS_ref will not be used and is only a dummy

!-----------------------------------------------------------------------------------
    
    ! Settings for the "dynamic" features
    
    ! Use all possible cluster-aerosol interactions
    ! If this setting conflicts with your cluster equations files, something should crash - but it's still good to not use wrong settings
    l_dynamic = .true.
    
    ! Simplification:
    ! Return the scavenged concentration for each vapor species onto each aerosol bin;
    ! here, aerosol growth by scavenging should be treated similarly to condensation
    coag_test_name = 'coag_molec'
    
    ! Accurate output:
    ! Return the scavenged concentration for each cluster onto each aerosol bin,
    ! as well as the cluster compositions;
    ! here, aerosol growth by scavenging should be treated similarly to coagulation
    !coag_test_name = 'coag_clust'
    
!-----------------------------------------------------------------------------------
!------------------------ End of user-defined input section ------------------------
!-----------------------------------------------------------------------------------
    
    ! A helper variable to make some output neater
    allocate(str_vapor(n_clustering_vapors))
    do i = 1,n_clustering_vapors
        str_vapor(i)(:) = '['//trim(names_vapor(i))//'] (cm^-3)'
    end do
    
    ! Allocate and initialize parameters for the dynamic set-up
    
    if (l_dynamic) then
        
        ! Coagulation from the cluster regime
        if (trim(coag_test_name) .eq. 'coag_molec') then
            
            allocate(c_coag_aero(naero,n_clustering_vapors))
        
        else if (trim(coag_test_name) .eq. 'coag_clust') then
                
            call get_system_size(nclust_syst=nclust_syst)
            
            allocate(c_coag_aero(naero,nclust_syst))
            allocate(clust_molec(nclust_syst,n_clustering_vapors))
            
        end if
        
        ! Size-classified flux of new particles out of the cluster regime
        allocate(c_out_bin(naero))
        allocate(comp_out_bin(naero,n_clustering_vapors))
        
        ! Composition of aerosols that shrink to the cluster regime
        allocate(nmols_evap(n_clustering_vapors))
        
    end if
    
!-----------------------------------------------------------------------------------
    
    ! Run a demo simulation for given time step(s)
    
    time_res = 1.d0*60.d0 ! s
    t = 0.d0
    
    ! Cluster plugin: All in- and output in SI units
    
    do it = 1,2
        
        ! AEROSOL MODEL I/O: The host model gives c_vapor
        
        if (l_dynamic) then
            
            ! Insert possible evaporating aerosol particles that shrink back to the cluster regime
            
            ! AEROSOL MODEL I/O: The host model gives c_evap and nmols_evap; NOTE: species not included in the cluster system must be transferred to gas phase by the aerosol model
            ! Demonstration: assume that the number and composition of evaporating aerosols is as below
            c_evap = 1.d5*1.d6
            if (it .eq. 1) c_evap = 0.d0
            nmols_evap = (/5.d0, 5.d0/)
            
            ! Aerosol parameters for the cluster-aerosol dynamics coupling
            
            ! NOTE: Every time that some of these parameters change, they must be updated by including them in the cluster dynamics routine call;
            ! only the updated parameters need to be inputted (but of course it's no harm to input also the unchanged parameters)
            
            ! Demonstration: assume that the size distribution of aerosols is as below
            
            ! Get bin diameter limits and mean diameters for the given number of bins (in practice this comes from the aerosol model)
            ! AEROSOL MODEL I/O: The host model gives naero, dp_aero_lim, dp_aero
            call dp_lim_aero(naero,dp_aero_lim,dp_aero)
            
            ! Get the aerosol number concentrations (in practice this comes from the aerosol model)
            ! AEROSOL MODEL I/O: The host model gives c_aero
            call lognormal_distr(ctot_aero,dmed_aero,sigma_aero,naero,c_aero,dp_aero)
            
            if (trim(coag_test_name) .eq. 'coag_molec') then
                
                call cluster_dynamics(names_vapor,c_vapor,cs_ref,temp,ipr,time_res,0.d0,&
                &    j_acdc,diameter_acdc,c_out=c_out,&
                &    naero=naero,dp_aero_lim=dp_aero_lim,dp_aero=dp_aero,c_aero=c_aero,&
                &    c_evap=c_evap,nmols_evap=nmols_evap,&
                &    c_coag_molec=c_coag_aero,&
                &    c_out_bin=c_out_bin,comp_out_bin=comp_out_bin)
                
            else if (trim(coag_test_name) .eq. 'coag_clust') then
                
                call cluster_dynamics(names_vapor,c_vapor,cs_ref,temp,ipr,time_res,0.d0,&
                &    j_acdc,diameter_acdc,c_out=c_out,&
                &    naero=naero,dp_aero_lim=dp_aero_lim,dp_aero=dp_aero,c_aero=c_aero,&
                &    c_evap=c_evap,nmols_evap=nmols_evap,&
                &    c_coag_clust=c_coag_aero,clust_molec=clust_molec,&
                &    c_out_bin=c_out_bin,comp_out_bin=comp_out_bin)
                
            end if
            
            ! Print some test output
            
            write(*,*)
            write(*,*) '---------- ClusterIn output after the given time step: ----------'
            write(*,*)
            write(*,*) 'Vapor concentrations'
            write(*,'(*(a15,1x))') str_vapor
            write(*,fmt_es) c_vapor*1.d-6
            write(*,*)
            write(*,*) 'Total concentration of new particles grown out of the cluster regime'
            write(char_tmp,fmt_es) j_acdc*1.d-6
            write(*,*) '(Average formation rate ', trim(char_tmp), ' cm^-3 s^-1)'
            write(char_tmp,fmt_es) c_out*1.d-6
            write(*,*) trim(char_tmp), ' cm^-3'
            write(*,*)
            write(*,*) 'Concentration grown to each aerosol bin and the average composition'
            do n = 1,naero
                if (c_out_bin(n) .gt. 0.d0) then
                    write(char_tmp,fmt_es) c_out_bin(n)*1.d-6
                    write(*,*) 'Bin ', n, ': ', trim(char_tmp), ' cm^-3, ', comp_out_bin(n,:), ' molecules'
                end if
            end do
            write(*,*)
            write(*,*) 'Concentration transferred to aerosol bins through cluster scavenging'
            do n = 1,naero
                if (any(c_coag_aero(n,:) .gt. 0.d0)) then
                    if (trim(coag_test_name) .eq. 'coag_molec') then
                        write(char_tmp,fmt_es) c_coag_aero(n,:)*1.d-6
                        write(*,*) 'Bin ', n, ': ', trim(char_tmp), ' vapor molecules (each species) cm^-3'
                    else if (trim(coag_test_name) .eq. 'coag_clust') then
                        write(char_tmp,fmt_es) sum(c_coag_aero(n,:))*1.d-6
                        write(*,*) 'Bin ', n, ': ', trim(char_tmp), ' clusters cm^-3'
                    end if
                end if
            end do
            
        else
            
            ! "Traditional" implementation in which only the formation rate is passed to the aerosol model
            
            call cluster_dynamics(names_vapor,c_vapor,cs_ref,temp,ipr,time_res,0.d0,&
            &    j_acdc,diameter_acdc,c_out=c_out)
            
            write(*,*)
            write(*,*) '---------- Formation rate output after the given time step: ----------'
            write(*,*)
            write(*,'(*(a20,1x))') str_vapor, 'CS_ref (s^-1)', 'T (K)', 'IPR (cm^-3 s^-1)', 'J (cm^-3 s^-1)'
            write(*,'(*(es20.5e4,1x))') c_vapor*1.d-6, cs_ref, temp, ipr*1.d-6, j_acdc*1.d-6
            
        end if
        
        ! Update the vapor concentrations in the aerosol model after the clustering routine
        ! AEROSOL MODEL I/O: The host model takes in updated c_vapor
        
        ! Finally, include the concentrations and compositions of formed and/or scavenged clusters wherever your aerosol model needs to have them:
        ! AEROSOL MODEL I/O: The host model updates the size distribution and composition
        ! Insert scavenged material to aerosol bins:
        ! c_coag_aero, or c_coag_aero and clust_molec
        ! Insert new particles and their composition to correct aerosol bin(s):
        ! c_out_bin
        ! comp_out_bin
        
        t = t + time_res
        
    end do
    
end program run_clusterin_example
