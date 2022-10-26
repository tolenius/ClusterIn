program run_clusterin_example_nchem

use acdc_aerosol_parameters, only : ctot_aero, dmed_aero, sigma_aero    ! Aerosol parameters for testing the dynamic coupling
use acdc_aerosol_parameters, only : dp_lim_aero, lognormal_distr        ! (will in practice be given by the aerosol dynamics model)
use acdc_aerosol_parameters, only : naero => naero_fixed

! use clusterin, only : cluster_dynamics                                 ! The main module for cluster input / output
! use acdc_simulation_setup, only : get_system_size                       ! For testing the dynamic features
include 'cluster_chem_use.inc'                                          ! For separate numbered cluster systems, one can get the modules from an .inc file

    implicit none
    
    ! All variables in SI units
    
    character(len=11), dimension(:,:), allocatable :: names_vapor_syst ! Vapor names for clustering species
    real(kind(1.d0)), allocatable :: c_vapor_tmp(:)                 ! Vapor concentrations for clustering species (m^-3)
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
    integer, allocatable :: nclust_syst(:)                          ! Number of clusters, possibly needed for testing the dynamic features
    real(kind(1.d0)) :: dp_aero_lim(naero+1), dp_aero(naero), c_aero(naero)  ! Aerosol parameters, needed for the dynamic features
    real(kind(1.d0)), allocatable  :: c_coag_aero_tmp(:,:)          ! Concentrations of coagulated molecules or clusters on aerosol bins (1/m^3)
    integer, allocatable  :: clust_molec_tmp(:,:)                   ! Cluster composition as numbers of vapor molecules (not specifying their possible charge)
    
    ! Formed particles that may be of different sizes (due to formation through different cluster-molecule and cluster-cluster collisions)
    real(kind(1.d0)), allocatable  :: c_out_bin(:)                  ! Concentrations of outgrown clusters to different aerosol bins (1/m^3)
    real(kind(1.d0)), allocatable :: comp_out_bin_tmp(:,:)          ! Composition of outgrown clusters to different aerosol bins (molec; note: real, not integer)
    
    ! Evaporation of larger particles back to the cluster regime
    real(kind(1.d0)) :: c_evap, c_evap_tmp                          ! Concentration of particles evaporating back to cluster regime (m^-3)
    real(kind(1.d0)), allocatable :: nmols_evap_tmp(:)              ! Composition of evaporating particles (molec; note: real, not integer)
    
    ! Misc.
    
    real(kind(1.d0)) :: time_res, t                                 ! Simulation times (s)
    
    ! Example of a host aerosol model system with n_spec_example vapor species that may or may not participate in clustering
    integer, parameter :: n_spec_example = 4
    character(len=11), dimension(n_spec_example) :: names_vapor_all
    real(kind(1.d0)) :: c_vapor_all(n_spec_example), nmols_evap_all(n_spec_example)
    logical :: l_clustering(n_spec_example) = .false.
    
    ! Numbers and indices of clustering systems
    integer :: n_clustering_syst, n_vapor_max                       ! Number of separate cluster systems and maximum number of vapors in any system
    integer, allocatable :: n_vapor_syst(:), ind_vapor_syst(:,:)    ! Number of vapors in each system and the indices of the vapors in c_vapor_all
    integer :: nclust_syst_max                                      ! Maximum number of clusters in any system, possibly needed for testing the dynamic features
    
    ! Some helper variables
    logical :: l_dynamic, l_clustering_numbers, l_found
    real(kind(1.d0)), allocatable  :: arr_tmp(:)
    character(len=11), dimension(:), allocatable :: str_vapor
    character(len=20) :: fmt_es = '(*(es15.2e4,1x))'
    character(len=200) :: char_tmp, fn_tmp, dir_clust_spec, fpre_clust_spec
    integer :: i, j, k, n, it, ind_v, ind_c, k_evap
    
    
    ! Input
    
    ! Simulation conditions
    
    ! Demonstration for how to get correct vapors from an array; particularly useful when using several cluster systems
    ! Let's say that these are the available vapors that may be included in the cluster systems
    names_vapor_all(1)(:) = 'O'
    names_vapor_all(2)(:) = 'A'
    names_vapor_all(3)(:) = 'N'
    names_vapor_all(4)(:) = 'D'
    c_vapor_all = (/1.d10, 2.d6, 1.d9, 1.d6/)*1.d6
    
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
    
    ! Clustering vapors
    
    ! Directory and file containing the clustering vapor names (one for each system, if there are several); used to automatically determine
    ! the vapors that are included in the cluster dynamics routines
    
    ! The file is assumed to be dir_clust_spec/fpre_clust_spec.txt, or dir_clust_spec_1/fpre_clust_spec_1.txt etc.
    !dir_clust_spec='.'
    !fpre_clust_spec='cluster_chem_spec_tmp'
    ! Separate numbered cluster systems
    dir_clust_spec='cluster_chem'
    fpre_clust_spec='cluster_chem_spec'
    
!-----------------------------------------------------------------------------------
!------------------------ End of user-defined input section ------------------------
!-----------------------------------------------------------------------------------

    ! Find information on the systems
    
    ! Find how many cluster systems there are, and how many vapor species they may include
    k = 0
    n_vapor_max = 0
    l_clustering_numbers = .true.
    l_found = .true.
    
    do while (l_found .and. l_clustering_numbers)
        
        k = k+1
        
        ! If there's only one system, it may not have any number suffix
        fn_tmp = trim(dir_clust_spec)//'/'//trim(fpre_clust_spec)//'.txt'
        inquire(file=fn_tmp,exist=l_found)
        
        if (l_found) then
            l_clustering_numbers = .false.
        else
            ! If no file is found, try with number suffices; note that one could also try to search e.g. for compiled .o files to ensure that
            ! the systems that are found are actually included in the code
            write(char_tmp,'(A,I0)') '_',k
            fn_tmp = trim(dir_clust_spec)//trim(char_tmp)//'/'//trim(fpre_clust_spec)//trim(char_tmp)//'.txt'
            inquire(file=fn_tmp,exist=l_found)
        end if
        
        if (l_found) then
            n_clustering_syst = k
            open(10,file=fn_tmp,action='read')
            ! Get the number of vapor species in the system from line 1
            read (10,*) i
            n_vapor_max = max(n_vapor_max,i)
            close(10)
        end if
        
    end do
    
    if (n_clustering_syst .ge. 1) then
        write(*,*)
        write(*,'(A,I0,A)') 'Found ',n_clustering_syst,' cluster dynamics system(s)'
        write(*,*)
    else
        write(*,*) 'Could not find any cluster system directories or files'
        stop
    end if
    
    ! Info on the cluster systems and their dimensions
    allocate(names_vapor_syst(n_vapor_max,n_clustering_syst))
    allocate(n_vapor_syst(n_clustering_syst))
    allocate(ind_vapor_syst(n_vapor_max,n_clustering_syst))
    
    ! Temporary helper variables that are used later
    allocate(c_vapor_tmp(n_vapor_max))
    allocate(str_vapor(n_vapor_max))
    allocate(arr_tmp(n_clustering_syst))
    
    ind_vapor_syst = 0
    
    do k = 1,n_clustering_syst
        
        ! Read the vapor names
        if (l_clustering_numbers) then
            write(char_tmp,'(A,I0)') '_',k
            fn_tmp = trim(dir_clust_spec)//trim(char_tmp)//'/'//trim(fpre_clust_spec)//trim(char_tmp)//'.txt'
        else
            fn_tmp = trim(dir_clust_spec)//'/'//trim(fpre_clust_spec)//'.txt'
        end if
        
        open(10,file=fn_tmp,action='read')
        ! Get the number of vapor species and the vapor names
        read (10,*) n_vapor_syst(k)
        i = 0
        do while (i .lt. n_vapor_syst(k))
            i = i+1
            read (10,*) names_vapor_syst(i,k)
        end do
        close(10)
        
        ! Find the clustering vapor indices in the array that contains all vapor names
        do i = 1,n_vapor_syst(k)
            
            l_found = .false.
            
            do j = 1,size(names_vapor_all)
                if (trim(names_vapor_syst(i,k)) .eq. trim(names_vapor_all(j))) then
                    ind_vapor_syst(i,k) = j
                    l_found = .true.
                    l_clustering(j) = .true.
                    !write(*,*) names_vapor_syst(i,k), ': species no. ', ind_vapor_syst(i,k)
                    exit
                end if
            end do
            
            if (.not. l_found) then
                write(*,*) 'Clustering vapor species ', names_vapor_syst(i,k), ' is not included in the host model'
                stop
            end if
            
        end do
        
    end do
    
!-----------------------------------------------------------------------------------
    
    ! Allocate and initialize parameters for the dynamic set-up
    
    if (l_dynamic) then
        
        ! Coagulation from the cluster regime
        if (trim(coag_test_name) .eq. 'coag_molec') then
            
            allocate(c_coag_aero_tmp(naero,n_vapor_max))
        
        else if (trim(coag_test_name) .eq. 'coag_clust') then
            
            allocate(nclust_syst(n_clustering_syst))
            
            do k = 1,n_clustering_syst
                
                ! call get_system_size(nclust_syst=nclust_syst(k))
                
                ! AEROSOL MODEL I/O: The example below is for including systems 1 and 2; modify if needed
                
                ! Separate numbered cluster systems
                if (k .eq. 1) then
                    call get_system_size_1(nclust_syst=nclust_syst(k))
                else if (k .eq. 2) then
                    call get_system_size_2(nclust_syst=nclust_syst(k))
                end if
                
            end do
            
            ! Maximum number of clusters in any system
            nclust_syst_max = maxval(nclust_syst)
            
            allocate(c_coag_aero_tmp(naero,nclust_syst_max))
            allocate(clust_molec_tmp(nclust_syst_max,n_vapor_max))
            
        end if
        
        ! Size-classified flux of new particles out of the cluster regime
        allocate(c_out_bin(naero))
        allocate(comp_out_bin_tmp(naero,n_vapor_max))
        
        ! Composition of aerosols that shrink to the cluster regime
        allocate(nmols_evap_tmp(n_vapor_max))
        
    end if
    
!-----------------------------------------------------------------------------------
    
    ! Run a demo simulation for given time step(s)
    
    time_res = 1.d0*60.d0 ! s
    t = 0.d0
    
    ! Cluster plugin: All in- and output in SI units
    
    do it = 1,2
        
        if (l_dynamic) then
            
            ! Insert possible evaporating aerosol particles that shrink back to the cluster regime
            
            ! AEROSOL MODEL I/O: The host model gives c_evap and nmols_evap_all
            ! Demonstration: assume that the number and composition of evaporating aerosols is as below
            c_evap = 1.d5*1.d6
            if (it .eq. 1) c_evap = 0.d0
            nmols_evap_all = (/2.d0, 5.d0, 5.d0, 1.d0/)
            
            ! Transfer to gas phase the species that are not included in any cluster system
            ! AEROSOL MODEL I/O: The host model takes in updated c_vapor_all and nmols_evap_all
            where (.not. l_clustering) c_vapor_all = c_vapor_all + c_evap*nmols_evap_all
            where (.not. l_clustering) nmols_evap_all = 0.d0
            
            ! The remaining core is returned to the cluster regime by inserting it in a suitable cluster system
            ! Index of the system
            k_evap = 0
            
            if (c_evap .gt. 0.d0 .and. maxval(nmols_evap_all) .gt. 0.d0) then
                
                do k = 1,n_clustering_syst
                
                    l_found = .true.
                    do i = 1,size(nmols_evap_all)
                        ! Test if the evaporating aerosol contains species that are not included in this system
                        if ((nmols_evap_all(i) .gt. 0.d0) .and. (.not. any(ind_vapor_syst(1:n_vapor_syst(k),k) .eq. i))) then
                            l_found = .false.
                            exit
                        end if
                    end do
                    
                    ! Select the first system that includes all the species in the evaporating aerosol
                    if (l_found) then
                        k_evap = k
                        exit
                    end if
                    
                end do
                
                ! If there is no system that contains all the species in the evaporating aerosol, select the system whose
                ! species have the highest molecular composition in the aerosol, and transfer other compounds to gas phase
                if (k_evap .eq. 0) then
                    
                    arr_tmp = 0.d0
                    do k = 1,n_clustering_syst
                        arr_tmp(k) = sum(nmols_evap_all(ind_vapor_syst(1:n_vapor_syst(k),k)))
                    end do
                    
                    k_evap = maxloc(arr_tmp,1)
                    
                    do i = 1,size(nmols_evap_all)
                        if ((nmols_evap_all(i) .gt. 0.d0) .and. &
                        &    (.not. any(ind_vapor_syst(1:n_vapor_syst(k_evap),k_evap) .eq. i))) then
                            !write(*,*) 'Removing species ',names_vapor_all(i)
                            ! AEROSOL MODEL I/O: The host model takes in updated c_vapor_all and nmols_evap_all
                            c_vapor_all(i) = c_vapor_all(i) + c_evap*nmols_evap_all(i)
                            nmols_evap_all(i) = 0.d0
                        end if
                    end do
                    
                end if
                
                write(*,*)
                write(*,*) 'Evaporating aerosol after processing:'
                write(*,*) 'The following composition will be inserted in cluster system no. ',k_evap,&
                &    ' (species ',names_vapor_syst(1:n_vapor_syst(k_evap),k_evap),')'
                write(*,*) names_vapor_all
                write(*,*) nmols_evap_all
                write(*,*)
            
            end if
            
        end if
        
        ! Simulate each clustering system, if there are several
        
        do k = 1,n_clustering_syst
            
            ! Number of vapors in this system; the temporary variables will be filled up to this index
            ind_v = n_vapor_syst(k)
            
            ! Include the clustering vapors that belong to this system
            ! AEROSOL MODEL I/O: The host model gives c_vapor_all
            c_vapor_tmp = 0.d0
            c_vapor_tmp(1:ind_v) = c_vapor_all(ind_vapor_syst(1:ind_v,k))
            
            ! A helper variable to make some output neater
            do i = 1,ind_v
                str_vapor(i)(:) = '['//trim(names_vapor_syst(i,k))//'] (cm^-3)'
            end do
            
            if (l_dynamic) then
                
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
                
                c_coag_aero_tmp = 0.d0
                comp_out_bin_tmp = 0.d0
                
                ! If this is the system chosen for shrinking aerosols, insert them here
                c_evap_tmp = 0.d0
                nmols_evap_tmp = 0.d0
                if (k .eq. k_evap) then
                    c_evap_tmp = c_evap
                    nmols_evap_tmp(1:ind_v) = nmols_evap_all(ind_vapor_syst(1:ind_v,k))
                end if
                
                if (trim(coag_test_name) .eq. 'coag_molec') then
                    
                    ! call cluster_dynamics(names_vapor_syst(1:ind_v,k),c_vapor_tmp(1:ind_v),cs_ref,temp,ipr,time_res,0.d0,&
                    ! &    j_acdc,diameter_acdc,c_out=c_out,&
                    ! &    naero=naero,dp_aero_lim=dp_aero_lim,dp_aero=dp_aero,c_aero=c_aero,&
                    ! &    c_evap=c_evap_tmp,nmols_evap=nmols_evap_tmp(1:ind_v),&
                    ! &    c_coag_molec=c_coag_aero_tmp(:,1:ind_v),&
                    ! &    c_out_bin=c_out_bin,comp_out_bin=comp_out_bin_tmp(:,1:ind_v))
                    
                    ! Choose which cluster system to call based on the system number - this is a bit clumsy, but since
                    ! routine names can't be parsed in Fortran it's the most straight-forward approach, for now
                    
                    ! AEROSOL MODEL I/O: The example below is for including systems 1 and 2; modify if needed
                    
                    ! Separate numbered cluster systems
                    if (k .eq. 1) then
                    
                        call cluster_dynamics_1(names_vapor_syst(1:ind_v,k),c_vapor_tmp(1:ind_v),cs_ref,temp,ipr,time_res,0.d0,&
                        &    j_acdc,diameter_acdc,c_out=c_out,&
                        &    naero=naero,dp_aero_lim=dp_aero_lim,dp_aero=dp_aero,c_aero=c_aero,&
                        &    c_evap=c_evap_tmp,nmols_evap=nmols_evap_tmp(1:ind_v),&
                        &    c_coag_molec=c_coag_aero_tmp(:,1:ind_v),&
                        &    c_out_bin=c_out_bin,comp_out_bin=comp_out_bin_tmp(:,1:ind_v))
                    
                    else if (k .eq. 2) then
                        
                        call cluster_dynamics_2(names_vapor_syst(1:ind_v,k),c_vapor_tmp(1:ind_v),cs_ref,temp,ipr,time_res,0.d0,&
                        &    j_acdc,diameter_acdc,c_out=c_out,&
                        &    naero=naero,dp_aero_lim=dp_aero_lim,dp_aero=dp_aero,c_aero=c_aero,&
                        &    c_evap=c_evap_tmp,nmols_evap=nmols_evap_tmp(1:ind_v),&
                        &    c_coag_molec=c_coag_aero_tmp(:,1:ind_v),&
                        &    c_out_bin=c_out_bin,comp_out_bin=comp_out_bin_tmp(:,1:ind_v))
                        
                    end if
                    
                else if (trim(coag_test_name) .eq. 'coag_clust') then
                    
                    ! Number of clusters in this system; the temporary variables will be filled up to this index
                    ind_c = nclust_syst(k)
                    
                    clust_molec_tmp = 0
                    
                    ! call cluster_dynamics(names_vapor_syst(1:ind_v,k),c_vapor_tmp(1:ind_v),cs_ref,temp,ipr,time_res,0.d0,&
                    ! &    j_acdc,diameter_acdc,c_out=c_out,&
                    ! &    naero=naero,dp_aero_lim=dp_aero_lim,dp_aero=dp_aero,c_aero=c_aero,&
                    ! &    c_evap=c_evap_tmp,nmols_evap=nmols_evap_tmp(1:ind_v),&
                    ! &    c_coag_clust=c_coag_aero_tmp(:,1:ind_c),clust_molec=clust_molec_tmp(1:ind_c,1:ind_v),&
                    ! &    c_out_bin=c_out_bin,comp_out_bin=comp_out_bin_tmp(:,1:ind_v))
                    
                    ! AEROSOL MODEL I/O: The example below is for including systems 1 and 2; modify if needed
                    
                    ! Separate numbered cluster systems
                    if (k .eq. 1) then
                        
                        call cluster_dynamics_1(names_vapor_syst(1:ind_v,k),c_vapor_tmp(1:ind_v),cs_ref,temp,ipr,time_res,0.d0,&
                        &    j_acdc,diameter_acdc,c_out=c_out,&
                        &    naero=naero,dp_aero_lim=dp_aero_lim,dp_aero=dp_aero,c_aero=c_aero,&
                        &    c_evap=c_evap_tmp,nmols_evap=nmols_evap_tmp(1:ind_v),&
                        &    c_coag_clust=c_coag_aero_tmp(:,1:ind_c),clust_molec=clust_molec_tmp(1:ind_c,1:ind_v),&
                        &    c_out_bin=c_out_bin,comp_out_bin=comp_out_bin_tmp(:,1:ind_v))
                        
                    else if (k .eq. 2) then
                        
                        call cluster_dynamics_2(names_vapor_syst(1:ind_v,k),c_vapor_tmp(1:ind_v),cs_ref,temp,ipr,time_res,0.d0,&
                        &    j_acdc,diameter_acdc,c_out=c_out,&
                        &    naero=naero,dp_aero_lim=dp_aero_lim,dp_aero=dp_aero,c_aero=c_aero,&
                        &    c_evap=c_evap_tmp,nmols_evap=nmols_evap_tmp(1:ind_v),&
                        &    c_coag_clust=c_coag_aero_tmp(:,1:ind_c),clust_molec=clust_molec_tmp(1:ind_c,1:ind_v),&
                        &    c_out_bin=c_out_bin,comp_out_bin=comp_out_bin_tmp(:,1:ind_v))
                        
                    end if
                    
                end if
                
                ! Print some test output
                
                write(*,*)
                write(*,*) '---------- ClusterIn output after the given time step: ----------'
                if (n_clustering_syst .gt. 1) then
                    write(*,'(A,I0,A)') '------------------- Cluster system no. ',k,' -------------------'
                end if
                write(*,*)
                write(*,*) 'Vapor concentrations'
                write(*,'(*(a15,1x))') str_vapor(1:ind_v)
                write(*,fmt_es) c_vapor_tmp(1:ind_v)*1.d-6
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
                        write(*,*) 'Bin ', n, ': ', trim(char_tmp), ' cm^-3, ', comp_out_bin_tmp(n,1:ind_v), ' molecules'
                    end if
                end do
                write(*,*)
                write(*,*) 'Concentration transferred to aerosol bins through cluster scavenging'
                do n = 1,naero
                    if (any(c_coag_aero_tmp(n,:) .gt. 0.d0)) then
                        if (trim(coag_test_name) .eq. 'coag_molec') then
                            write(char_tmp,fmt_es) c_coag_aero_tmp(n,1:ind_v)*1.d-6
                            write(*,*) 'Bin ', n, ': ', trim(char_tmp), ' vapor molecules (each species) cm^-3'
                        else if (trim(coag_test_name) .eq. 'coag_clust') then
                            write(char_tmp,fmt_es) sum(c_coag_aero_tmp(n,1:ind_c))*1.d-6
                            write(*,*) 'Bin ', n, ': ', trim(char_tmp), ' clusters cm^-3'
                        end if
                    end if
                end do
                
            else
                
                ! "Traditional" implementation in which only the formation rate is passed to the aerosol model
                
                ! call cluster_dynamics(names_vapor_syst(1:ind_v,k),c_vapor_tmp(1:ind_v),cs_ref,temp,ipr,time_res,0.d0,&
                ! &    j_acdc,diameter_acdc,c_out=c_out)
                
                ! AEROSOL MODEL I/O: The example below is for including systems 1 and 2; modify if needed
                
                ! Separate numbered cluster systems
                if (k .eq. 1) then
                    
                    call cluster_dynamics_1(names_vapor_syst(1:ind_v,k),c_vapor_tmp(1:ind_v),cs_ref,temp,ipr,time_res,0.d0,&
                    &    j_acdc,diameter_acdc,c_out=c_out)
                    
                else if (k .eq. 2) then
                    
                    call cluster_dynamics_2(names_vapor_syst(1:ind_v,k),c_vapor_tmp(1:ind_v),cs_ref,temp,ipr,time_res,0.d0,&
                    &    j_acdc,diameter_acdc,c_out=c_out)
                    
                end if
                
                write(*,*)
                write(*,*) '---------- Formation rate output after the given time step: ----------'
                if (n_clustering_syst .gt. 1) then
                    write(*,'(A,I0,A)') '------------------ Cluster system no. ',k,' ------------------'
                end if
                write(*,*)
                write(*,'(*(a20,1x))') str_vapor(1:ind_v), 'CS_ref (s^-1)', 'T (K)', 'IPR (cm^-3 s^-1)', 'J (cm^-3 s^-1)'
                write(*,'(*(es20.5e4,1x))') c_vapor_tmp(1:ind_v)*1.d-6, cs_ref, temp, ipr*1.d-6, j_acdc*1.d-6
                
            end if
            
            ! Update the vapor concentrations in the aerosol model after the clustering routine
            ! AEROSOL MODEL I/O: The host model takes in updated c_vapor_all
            c_vapor_all(ind_vapor_syst(1:ind_v,k)) = c_vapor_tmp(1:ind_v)
            
            ! Finally, include the concentrations and compositions of formed and/or scavenged clusters wherever your aerosol model needs to have them:
            ! AEROSOL MODEL I/O: The host model updates the size distribution and composition
            ! Insert scavenged material to aerosol bins:
            ! c_coag_aero_tmp(:,1:ind_v), or c_coag_aero_tmp(:,1:ind_c) and clust_molec_tmp(1:ind_c,1:ind_v)
            ! Insert new particles and their composition to correct aerosol bin(s):
            ! c_out_bin
            ! comp_out_bin_tmp(:,1:ind_v)
            
        end do
        
        t = t + time_res
        
    end do
    
end program run_clusterin_example_nchem
