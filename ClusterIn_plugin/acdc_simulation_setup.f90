module acdc_simulation_setup
use acdc_system

implicit none

logical, parameter :: solve_ss = .false.        ! Solve the steady state or run only for a given time
logical, parameter :: j_flux = .false.          ! J calculated as a flux or as particles per time step; note that the former
                                                ! is automatically used for solve_ss = .true. regardless of the j_flux logical

integer, parameter :: nbins = 5                 ! Number of size bins for size-classifying the clusters (if used)

real(kind(1.d0)), parameter :: pi=4.d0*atan(1.d0)


contains

!--------------------------------------------------------------------------------------------------
! Settings for the monomer concentrations
!--------------------------------------------------------------------------------------------------

subroutine sources_and_constants(neqn,source,isconst,fitted)

    implicit none
    
    real(kind(1.d0)) :: source(neqn)
    logical :: isconst(neqn)
    integer :: neqn, fitted(neqn,0:neqn), cluster_numbers(n_1A_clusters), n, i
    integer :: n_monomers(n_monomer_types)
    
    source = 0.d0                               ! Initialize all source terms to 0
    isconst = .false.                           ! Initialize all concentrations to vary freely according to the eqs.
    fitted = 0
    
    ! NOTE: If some concentrations are forced to be constant by excluding their equations (Perl option --no_eq),
    ! the below settings are NOT used


    ! Set some or all monomers to have constant or fitted concentrations
    
    ! If everything related to the constant or fitted concentrations is commented out,
    ! all concentrations will be explicitly determined from the equations
    
    if (solve_ss) then
    
        if (small_set_mode) then

            ! Use constant concentrations for electrically neutral monomers
            isconst(neutral_monomers) = .true.          ! Neutral monomer concentrations are constant
                                                        ! (Ionic monomers can vary according to the ion production rate)

            ! ! Fit concentrations (here for 1A)
            ! isconst(n1A) = .false.
            ! fitted(1,0) = 1                                        ! 1 concentration is fitted
            ! call clusters_with_1_A(cluster_numbers)
            ! n = n_1A_clusters
            ! fitted(2,0:n) = (/n1A, n-1, cluster_numbers(2:n)/)     ! Concentration of 1A is fitted using
                                                                   ! !  n-1 other concentrations: [1A]+... = const.
            ! ! Other monomer concentrations will remain as set above
        
        else
        
            ! Use constant concentrations for all monomers
            call monomer_indices(n_monomers)
            isconst(n_monomers) = .true.
        
        end if
    
    end if
    
end subroutine sources_and_constants

!--------------------------------------------------------------------------------------------------
! Parameters related to the size of the simulation system, including
! the numbers of clusters and equations, indices of outgoing fluxes, ...
!--------------------------------------------------------------------------------------------------

subroutine get_system_size(neq_syst,nclust_syst,nout_syst,nclust_out,nrange_coag,nrange_out,diameter_max_syst)
!use restriction_criteria, only : diameter_max_loop

    implicit none
    
    integer, intent(out), optional :: neq_syst, nclust_syst             ! Numbers of equations and clusters
    integer, intent(out), optional :: nout_syst(n_charges)              ! Indices of outgoing fluxes
    integer, intent(out), optional :: nclust_out                        ! Number of individual clusters growing out
    integer, intent(out), optional :: nrange_coag(2)                    ! Index range of scavenged clusters
    integer, intent(out), optional :: nrange_out(2)                     ! Index range of outgrown clusters
    real(kind(1.d0)), intent(out), optional :: diameter_max_syst        ! Max. diameter in system
    integer, save :: nclust_saved = 0
    integer, allocatable :: clust_comp(:,:)                             ! Cluster composition in the loop mode
    
    if (small_set_mode) then
        
        if (present(neq_syst)) neq_syst = neq
        if (present(nclust_syst)) nclust_syst = nclust
        if (present(nout_syst)) nout_syst = nout_all
        
        if (present(nclust_out)) nclust_out = nclust_out_saved
        if (present(nrange_coag)) nrange_coag = nrange_coag_saved
        if (present(nrange_out)) nrange_out = nrange_out_saved
        
        if (present(diameter_max_syst)) diameter_max_syst = diameter_max*1.d-9
        
    else
        
        if (n_mol_types .gt. 1) then
            allocate(clust_comp(nclust_max,n_mol_types))
            call get_molecule_numbers(nclust_saved,clust_comp)
        else
            nclust_saved = nclust_max
        end if
        
        if (present(neq_syst)) neq_syst = nclust_saved+(neq_max-nclust_max)
        if (present(nclust_syst)) nclust_syst = nclust_saved
        if (present(nout_syst)) nout_syst = nclust_saved+1    ! Loop mode: last, additional index
        
        !if (present(diameter_max_syst)) diameter_max_syst = diameter_max_loop
        
    end if
    
end subroutine get_system_size

!--------------------------------------------------------------------------------------------------
! Indices of the vapor molecules in the cluster and molecule type arrays
!--------------------------------------------------------------------------------------------------

subroutine get_vapor_indices(names_vapor,n1vapor,nmol_vapor)

    implicit none
    
    character(len=11), dimension(n_neutral_monomers), intent(in) :: names_vapor
    integer, intent(out) :: n1vapor(n_neutral_monomers), nmol_vapor(n_neutral_monomers)
    character(len=11), dimension(n_monomer_types) :: names_monomers
    integer :: n_monomers(n_monomer_types)
    character(len=11), dimension(n_mol_types) :: names_mol_types
    integer :: i, j
    logical :: lfound
    
    ! Find the cluster and molecule numbers for the included vapors
    
    n1vapor = 0
    nmol_vapor = 0
    
    call monomer_names(names_monomers)
    call monomer_indices(n_monomers)
    call molecule_names(names_mol_types)
    
    do i = 1,size(names_vapor)
    
        lfound = .false.
        do j = 1,size(names_monomers)
            if (trim('1'//names_vapor(i)(:)) .eq. trim(names_monomers(j)(:))) then
                n1vapor(i) = n_monomers(j)
                lfound = .true.
                !write(*,*) names_vapor(i)(:), ': cluster no. ', n1vapor(i)
                exit
            end if
        end do
        if (.not. lfound) then
            write(*,*) 'Monomer ', names_vapor(i)(:), ' not included in the ACDC set-up'
            stop
        end if
        
        lfound = .false.
        do j = 1,size(names_mol_types)
            if (trim(names_vapor(i)(:)) .eq. trim(names_mol_types(j)(:))) then
                nmol_vapor(i) = j
                lfound = .true.
                !write(*,*) names_vapor(i)(:), ': molecule no. ', nmol_vapor(i)
                exit
            end if
        end do
        if (.not. lfound) then
            write(*,*) 'Molecule ', names_vapor(i)(:), ' not included in the ACDC set-up'
            stop
        end if
        
    end do
    
end subroutine get_vapor_indices

!--------------------------------------------------------------------------------------------------
! Index of the boundary cluster (most) corresponding to the given composition, needed for
! including larger particles evaporating back to the cluster regime
!--------------------------------------------------------------------------------------------------

subroutine get_bound_index(nmol_vapor,molec_vapor,nclust_bound,diff_molec_vapor)

    implicit none
    
    integer, intent(in) :: nmol_vapor(n_neutral_monomers)                       ! Indices of the vapor molecules in the ACDC molecule type array
    real(kind(1.d0)), intent(in) :: molec_vapor(n_neutral_monomers)             ! Molecular composition of evaporating particles (note: real, not integer)
    integer, intent(out) :: nclust_bound                                        ! Index of the corresponding ACDC cluster
    real(kind(1.d0)), intent(out) :: diff_molec_vapor(n_neutral_monomers)       ! Difference in the molecular compositions of the particle and the ACDC cluster
    integer, save :: bound_clusters(n_bound), nmols_bound(n_bound,n_mol_types)
    real(kind(1.d0)), save :: frac_bound(n_bound,n_mol_types), diff_th
    real(kind(1.d0)) :: frac_vapor(n_neutral_monomers), frac_mol_types(n_mol_types), diff_sum(n_bound)
    integer, save :: ch_state_bound(n_bound)
    integer :: ch_state(nclust)
    integer :: i, ind(1)
    logical, save :: firstcall = .true.
    
    nclust_bound = 0
    
    if (firstcall) then
    
        firstcall = .false.
        
        call get_bound(bound_clusters,nmols_bound)
        
        ! Boundary cluster molar composition
        do i = 1,n_bound
            frac_bound(i,:) = real(nmols_bound(i,:),kind=kind(1.d0))/real(sum(nmols_bound(i,:)),kind=kind(1.d0))
        end do
        
        ! Boundary cluster charging state
        call get_charging_state(ch_state)
        ch_state_bound = ch_state(bound_clusters)
        
        ! Threshold for printing a warning about the difference between the given composition and
        ! the closest boundary composition, given as the deviation in the molar composition for each species
        diff_th = 5.d-2
        
    end if
    
    ! Incoming molar composition
    do i = 1,n_neutral_monomers
        frac_vapor(i) = molec_vapor(i)/sum(molec_vapor)
    end do
    
    !write(*,*) 'Aerosol molecular composition: ',molec_vapor
    !write(*,*) 'Aerosol molar composition: ',frac_vapor
    
    ! For simplicity, the aerosol dynamics model is assumed to give the vapor species, which are here mapped to
    ! correspond to the order of the ACDC molecule types (which may not be the same as monomer types)
    
    frac_mol_types = 0.d0
    
    do i = 1,size(nmol_vapor)
        frac_mol_types(nmol_vapor(i)) = frac_vapor(i)
    end do
    
    diff_sum = 0.d0
    
    do i = 1,n_bound
        diff_sum(i) = sum(abs(frac_bound(i,:)-frac_mol_types))
    end do
    
    ! For now, assume that the evaporating particles are electrically neutral
    ind = minloc(diff_sum, ch_state_bound .eq. 1)
    nclust_bound = bound_clusters(ind(1))
    
    ! Difference in the molecular composition
    diff_molec_vapor = molec_vapor - real(nmols_bound(ind(1),nmol_vapor),kind=kind(1.d0))
    
    !write(*,*) 'ACDC molecular composition: ',nmols_bound(ind(1),nmol_vapor),' (cluster ',nclust_bound,')'
    !write(*,*) 'ACDC molar composition: ',frac_bound(ind(1),nmol_vapor)
    !write(*,*) 'Molecular difference: ',diff_molec_vapor
    
    if (any(abs(frac_bound(ind(1),:)-frac_mol_types) .gt. diff_th)) then
        write(*,*) 'Warning: no ACDC boundary cluster corresponding closely to the evaporating particles'
        write(*,*) 'Evaporating particles'' molar composition ', frac_mol_types
        write(*,*) 'Closest boundary composition ', frac_bound(ind(1),:)
    end if

end subroutine get_bound_index

!--------------------------------------------------------------------------------------------------
! Scavenging sink for each cluster and scavenged molecules for each aerosol size bin based on
! the exact aerosol size distribution, needed for including coagulation from the cluster regime
! onto larger particles
!
! The following routines need to be provided by the aerosol model interface:
! get_cs_aero:          Subroutine that gives the sink coefficient of each cluster onto each aerosol bin
!
!--------------------------------------------------------------------------------------------------

subroutine get_scav_parameters(names_vapor,temperature,c_coag,cs_per_clust,&
    &    names_nocharge,molec_per_aero,clust_per_aero,nmols_nocharge_clust)

    use acdc_aerosol_parameters, only : get_cs_aero

    implicit none
    
    ! Optional input
    character(len=11), dimension(n_neutral_monomers), intent(in), optional :: names_vapor   ! Names of the vapor species, if the output is wanted in this order
    real(kind(1.d0)), intent(in), optional :: temperature                       ! Temperature (possibly needed for the sink rates) (K)
    real(kind(1.d0)), intent(inout), optional :: c_coag(:)                      ! Scavenged concentration for each cluster, to be mapped onto the aerosol bins (1/m^3)
    ! Optional output
    real(kind(1.d0)), intent(out), optional :: cs_per_clust(:)                  ! Total sink for each cluster; can be used within the ACDC simulation (1/s)
    character(len=11), dimension(n_neutral_monomers), intent(out), optional :: names_nocharge   ! Names of the molecules transfered by coagulation
    real(kind(1.d0)), intent(out), optional :: molec_per_aero(:,:)              ! Total coagulated molecules for each aerosol bin and vapor (molec/m^3);
                                                                                ! size: (naero,n_neutral_monomers)
    real(kind(1.d0)), intent(out), optional :: clust_per_aero(:,:)              ! Coagulated clusters for each aerosol bin and cluster composition (1/m^3);
                                                                                ! size: (naero,nclust)
    integer, intent(out), optional :: nmols_nocharge_clust(:,:)                 ! Cluster composition omitting charge; e.g. HSO4- counted as H2SO4;
                                                                                ! size: (nclust,n_neutral_monomers)
    
    integer, save :: nclust_syst = 0, naero = 0                                 ! Numbers of clusters and aerosol bins
    real(kind(1.d0)), allocatable, save :: dp_aero(:), mp_aero(:), c_aero(:)    ! Aerosol diameters (m), masses (kg) and concentrations (m^-3)
    real(kind(1.d0)), allocatable, save :: dp_clust(:), mp_clust(:)             ! Cluster diameters (m) and masses (kg)
    real(kind(1.d0)) :: pres                                                    ! Pressure (possibly needed for the sink rates) (Pa)
    real(kind(1.d0)), allocatable, save :: cs_per_clust_aero(:,:), cs_per_clust_saved(:)    ! Sink for each cluster onto each aerosol bin / total sink (1/s)
    character(len=11), dimension(n_neutral_monomers), save :: names_nocharge_acdc, names_nocharge_saved
    integer, allocatable, save :: nmols_nocharge_clust_acdc(:,:), nmols_nocharge_clust_saved(:,:), ind_nmols_nocharge_zero(:)
    integer, save :: naero_saved = 0, nclust_nmols_zero = 0
    integer, allocatable :: ind_tmp(:)
    integer :: i, j
    logical :: lfound
    logical, save :: firstcall = .true.
    
    if (firstcall) then
        
        firstcall = .false.
        
        call get_system_size(nclust_syst=nclust_syst)
        
        allocate(dp_clust(nclust_syst))
        allocate(mp_clust(nclust_syst))
        allocate(nmols_nocharge_clust_acdc(nclust_syst,n_neutral_monomers))
        allocate(nmols_nocharge_clust_saved(nclust_syst,n_neutral_monomers))
        
        allocate(cs_per_clust_saved(nclust_syst))
        
        ! Get cluster diameters and masses
        call get_diameter(dp_clust)
        call get_mass(mp_clust)
        dp_clust = dp_clust*1.d-9                   ! nm -> m
        mp_clust = mp_clust*1.d-3/6.02214179d23       ! amu -> kg
        
        ! Get cluster compositions in terms of summed neutral and charged molecules of each species
        call molecule_names_nocharge(names_nocharge_acdc)
        call get_nmols_nocharge(nmols_nocharge_clust_acdc)
        names_nocharge_saved = names_nocharge_acdc
        nmols_nocharge_clust_saved = nmols_nocharge_clust_acdc
        
        ! If some clusters (e.g. generic ions) contain no tracked vapor molecules, save their indices and
        ! set their coagulated concentration to zero, to be sure that they are not included anywhere
        allocate(ind_tmp(nclust_syst))
        ind_tmp = 0
        do i = 1,nclust_syst
            if (all(nmols_nocharge_clust_saved(i,:) .eq. 0)) then
                nclust_nmols_zero = nclust_nmols_zero+1
                ind_tmp(nclust_nmols_zero) = i
            end if
        end do
        if (nclust_nmols_zero .gt. 0) then
            allocate(ind_nmols_nocharge_zero(nclust_nmols_zero))
            ind_nmols_nocharge_zero = ind_tmp(1:nclust_nmols_zero)
        end if
        
    end if
    
    call get_ambient_params(naero=naero)
    
    if (naero .ne. naero_saved) then
        
        naero_saved = naero
        
        if (allocated(dp_aero)) deallocate(dp_aero)
        if (allocated(mp_aero)) deallocate(mp_aero)
        if (allocated(c_aero)) deallocate(c_aero)
        if (allocated(cs_per_clust_aero)) deallocate(cs_per_clust_aero)
        
        allocate(dp_aero(naero))
        allocate(mp_aero(naero))
        allocate(c_aero(naero))
        allocate(cs_per_clust_aero(naero,nclust_syst))
        
    end if
    
    ! Find the order of the non-charged molecules if the vapor names are given as input
    if (present(names_vapor)) then
        
        if (any(names_vapor .ne. names_nocharge_saved)) then
            
            nmols_nocharge_clust_saved = 0
            
            do i = 1,size(names_vapor)
            
                lfound = .false.
                do j = 1,size(names_nocharge_acdc)
                    if (trim(names_vapor(i)(:)) .eq. trim(names_nocharge_acdc(j)(:))) then
                        ! Re-order the species in the non-charged cluster composition array to obtain the output in the correct order
                        nmols_nocharge_clust_saved(:,i) = nmols_nocharge_clust_acdc(:,j)
                        lfound = .true.
                        write(*,*) names_vapor(i)(:), ': non-charged species no. ', j
                        exit
                    end if
                end do
                if (.not. lfound) then
                    write(*,*) 'Molecule ', names_vapor(i)(:), ' not found among the non-charged forms of the ACDC species'
                    stop
                end if
            
            end do
            
            names_nocharge_saved = names_vapor
            
        end if
        
    end if
    
    if (present(names_nocharge)) names_nocharge = names_nocharge_saved    
    if (present(nmols_nocharge_clust)) nmols_nocharge_clust = nmols_nocharge_clust_saved
    
    ! Get the sink coefficient for each cluster onto each aerosol bin according to a user-defined approach
    if (.not. present(temperature)) then
        write(*,*) 'get_scav_parameters: Temperature must be given as input at every call in the current implementation'
        stop
    end if
    call get_ambient_params(dp_aero=dp_aero,mp_aero=mp_aero,c_aero=c_aero,pres=pres)
    call get_cs_aero(nclust_syst,naero,dp_aero,mp_aero,c_aero,dp_clust,mp_clust,temperature,pres,cs_per_clust_aero)
    
    ! Calculate the total coagulation sink for each cluster
    cs_per_clust_saved = 0.d0
    do i = 1, nclust_syst
        cs_per_clust_saved(i) = sum(cs_per_clust_aero(:,i))
    end do
    
    if (present(cs_per_clust)) cs_per_clust = cs_per_clust_saved
    
    if (present(c_coag) .and. nclust_nmols_zero .gt. 0) then
        c_coag(ind_nmols_nocharge_zero) = 0.d0
    end if
    
    ! Calculate the total numbers of clusters coagulated onto each aerosol bin (1/m^3)
    if (present(clust_per_aero)) then
        
        clust_per_aero = 0.d0
        
        do i = 1, nclust_syst
            do j = 1, naero
                clust_per_aero(j,i) = cs_per_clust_aero(j,i)/cs_per_clust_saved(i)*c_coag(i)
            end do
        end do
        
    end if
    
    ! Calculate the total numbers of molecules of each vapor species coagulated onto each aerosol bin (molec/m^3)
    if (present(molec_per_aero)) then
        
        molec_per_aero = 0.d0
        
        do i = 1, nclust_syst
            do j = 1, naero
                molec_per_aero(j,:) = molec_per_aero(j,:) + cs_per_clust_aero(j,i)/cs_per_clust_saved(i)*&
                &    c_coag(i)*real(nmols_nocharge_clust_saved(i,:),kind=kind(1.d0))
            end do
        end do
        
    end if

end subroutine get_scav_parameters

!--------------------------------------------------------------------------------------------------
! Formed particles and their average composition per each aerosol size bin, based on the individual
! cluster compositions that grow out of the cluster regime, needed when the new particles that are
! inserted in the aerosol model may come in several different sizes
!--------------------------------------------------------------------------------------------------

subroutine get_jbin_parameters(c_out_clust,c_per_aero,comp_nocharge_per_aero,nmols_nocharge_out,&
    &    names_vapor,names_nocharge)
    
    implicit none
    
    ! Input
    real(kind(1.d0)), intent(in) :: c_out_clust(:)                              ! Concentration for each outgrown cluster, to be mapped onto the aerosol bins (1/m^3)
    ! Optional input
    character(len=11), dimension(n_neutral_monomers), intent(in), optional :: names_vapor   ! Names of the vapor species, if the output is wanted in this order
    ! Output
    real(kind(1.d0)), intent(out) :: c_per_aero(:)                              ! New particles for each aerosol bin (1/m^3);
                                                                                ! size: (naero)
    real(kind(1.d0)), intent(out) :: comp_nocharge_per_aero(:,:)                ! Average composition of new particles for each aerosol bin (molec);
                                                                                ! size: (naero,n_neutral_monomers)
    ! Optional output
    integer, intent(out), optional :: nmols_nocharge_out(:,:)                   ! Composition of outgrown clusters (molec);
                                                                                ! size: (nclust_out,n_neutral_monomers)
    character(len=11), dimension(n_neutral_monomers), intent(out), optional :: names_nocharge   ! Names of the molecules in the new particles
    
    integer, save :: nclust_out = 0, naero = 0                                  ! Numbers of outgrown clusters and aerosol bins
    real(kind(1.d0)), allocatable, save :: dp_clust_out(:)                      ! Diameters of outgrown clusters (m)
    integer, allocatable, save :: aero_for_clust_out(:)                         ! Aerosol bin index for each outgrown cluster
    real(kind(1.d0)), allocatable, save :: dp_aero_lim(:), dp_aero_lim_saved(:) ! Diameters corresponding to the aerosol bin limits (m)
    character(len=11), dimension(n_neutral_monomers), save :: names_nocharge_acdc, names_nocharge_saved
    integer, allocatable, save :: nmols_nocharge_out_acdc(:,:), nmols_nocharge_out_saved(:,:)   ! Outgrown cluster composition omitting charge; e.g. HSO4- counted as H2SO4
    integer, save :: naero_saved = 0
    integer :: i, j
    logical :: lfound
    logical, save :: firstcall = .true.
    
    if (firstcall) then
        
        firstcall = .false.
        
        call get_system_size(nclust_out=nclust_out)
        
        allocate(dp_clust_out(nclust_out))
        allocate(nmols_nocharge_out_acdc(nclust_out,n_neutral_monomers))
        allocate(nmols_nocharge_out_saved(nclust_out,n_neutral_monomers))
        
        allocate(aero_for_clust_out(nclust_out))
        
        ! Get cluster diameters
        call get_diameter_out(dp_clust_out)
        dp_clust_out = dp_clust_out*1.d-9           ! nm -> m
        
        ! Get cluster compositions in terms of summed neutral and charged molecules of each species
        call molecule_names_nocharge(names_nocharge_acdc)
        call get_nmols_nocharge_out(nmols_nocharge_out_acdc)
        names_nocharge_saved = names_nocharge_acdc
        nmols_nocharge_out_saved = nmols_nocharge_out_acdc
        
    end if
    
    call get_ambient_params(naero=naero)
    
    if (naero .ne. naero_saved) then
        
        naero_saved = naero
        
        if (allocated(dp_aero_lim)) deallocate(dp_aero_lim)
        allocate(dp_aero_lim(naero+1))
        
        if (allocated(dp_aero_lim_saved)) deallocate(dp_aero_lim_saved)
        allocate(dp_aero_lim_saved(naero+1))
        
        dp_aero_lim_saved = 0.d0
        
    end if
    
    ! Find the order of the non-charged molecules if the vapor names are given as input
    if (present(names_vapor)) then
        
        if (any(names_vapor .ne. names_nocharge_saved)) then
            
            nmols_nocharge_out_saved = 0
            
            do i = 1,size(names_vapor)
            
                lfound = .false.
                do j = 1,size(names_nocharge_acdc)
                    if (trim(names_vapor(i)(:)) .eq. trim(names_nocharge_acdc(j)(:))) then
                        ! Re-order the species in the non-charged cluster composition array to obtain the output in the correct order
                        nmols_nocharge_out_saved(:,i) = nmols_nocharge_out_acdc(:,j)
                        lfound = .true.
                        write(*,*) names_vapor(i)(:), ': non-charged species no. ', j
                        exit
                    end if
                end do
                if (.not. lfound) then
                    write(*,*) 'Molecule ', names_vapor(i)(:), ' not found among the non-charged forms of the ACDC species'
                    stop
                end if
            
            end do
            
            names_nocharge_saved = names_vapor
            
        end if
        
    end if
    
    if (present(names_nocharge)) names_nocharge = names_nocharge_saved
    
    call get_ambient_params(dp_aero_lim=dp_aero_lim)
    
    ! Update the cluster mapping if the aerosol bin limits have changed
    if (maxval(abs((dp_aero_lim-dp_aero_lim_saved)/max(dp_aero_lim_saved,1.d-20))) .gt. 1.d-3) then
        
        dp_aero_lim_saved = dp_aero_lim
        
        ! Map the outgrown clusters to the aerosol bins
        write(*,*)
        write(*,*) 'Determining aerosol bins for the formed new particles'
        write(*,*)
        
        aero_for_clust_out = 0
        
        do i = 1, nclust_out
            lfound = .false.
            do j = 1, naero
                if ((dp_clust_out(i) .ge. dp_aero_lim(j)) .and. (dp_clust_out(i) .lt. dp_aero_lim(j+1))) then
                    aero_for_clust_out(i) = j
                    lfound = .true.
                    !write(*,*) 'Outgrown cluster of diameter ', dp_clust_out(i)*1.d9, ' nm placed in aerosol bin no. ',j
                    exit
                end if
            end do
            if (.not. lfound) then
                write(*,*) 'Failed to determine the aerosol bin for a new particle of diameter ', dp_clust_out(i)*1.d9, ' nm'
                stop
            end if
        end do
        
        if (all(aero_for_clust_out .eq. 1)) then
            write(*,*) 'All clusters growing into aerosol bin no. 1'
            write(*,*) '- binning of J can still be good as it gives the exact average composition'
            write(*,*)
        else if (.not. any(aero_for_clust_out .eq. 1)) then
            write(*,*) 'Warning: No clusters growing into aerosol bin no. 1 - increase the size of the smallest bin?'
            write(*,*)
        end if
        
    end if
    
    ! Calculate the total numbers of new particles and their average composition for each aerosol bin (1/m^3)
    c_per_aero = 0.d0
    comp_nocharge_per_aero = 0.d0
    
    do i = 1, nclust_out
        j = aero_for_clust_out(i)
        c_per_aero(j) = c_per_aero(j) + c_out_clust(i)
        comp_nocharge_per_aero(j,:) = comp_nocharge_per_aero(j,:) + c_out_clust(i)*&
        &    real(nmols_nocharge_out_saved(i,:),kind=kind(1.d0))
    end do
    do j = 1, naero
        if (c_per_aero(j) .gt. 0.d0) comp_nocharge_per_aero(j,:) = comp_nocharge_per_aero(j,:)/c_per_aero(j)
    end do
    
    if (present(nmols_nocharge_out)) nmols_nocharge_out = nmols_nocharge_out_saved

end subroutine get_jbin_parameters

!--------------------------------------------------------------------------------------------------
! Saved aerosol model parameters (bin diameters, concentrations, ...), needed for routines that
! deal with scavenged or outgrown clusters - the routine only takes in and saves the parameters
! in order to have them accessible to all cluster dynamics programs
!--------------------------------------------------------------------------------------------------

subroutine get_ambient_params(naero_in,naero,dp_aero_lim_in,dp_aero_lim,dp_aero_in,dp_aero,&
    &    mp_aero_in,mp_aero,c_aero_in,c_aero,pres_in,pres)

    implicit none
    
    ! Aerosol parameters
    integer, intent(in), optional :: naero_in                           ! Number of aerosol size bins
    integer, intent(out), optional :: naero
    real(kind(1.d0)), intent(in), optional :: dp_aero_lim_in(:)         ! Bin limits (m)
    real(kind(1.d0)), intent(out), optional :: dp_aero_lim(:) 
    real(kind(1.d0)), intent(in), optional :: dp_aero_in(:)             ! Bin mean diameter (m)
    real(kind(1.d0)), intent(out), optional :: dp_aero(:)
    real(kind(1.d0)), intent(in), optional :: mp_aero_in(:)             ! Bin single-particle mass (kg)
    real(kind(1.d0)), intent(out), optional :: mp_aero(:)
    real(kind(1.d0)), intent(in), optional :: c_aero_in(:)              ! Size-binned aerosol concentration (m^-3)
    real(kind(1.d0)), intent(out), optional :: c_aero(:)
    
    ! Other ambient parameters
    real(kind(1.d0)), intent(in), optional :: pres_in                   ! Ambient pressure (Pa)
    real(kind(1.d0)), intent(out), optional :: pres
    
    integer, save :: naero_saved = 0
    real(kind(1.d0)), allocatable, save :: dp_aero_lim_saved(:), dp_aero_saved(:), mp_aero_saved(:), c_aero_saved(:)
    real(kind(1.d0)), save :: rhop_def = 1500.d0                        ! Default particle density, used to get the mass if no input is given (kg m^-3)
    real(kind(1.d0)), save :: pres_saved = 101325.d0                    ! Default pressure (Pa)
    
    integer, parameter :: n_params = 6                                  ! Number of storable parameters
    integer, save :: n_naero = 1, n_dp_aero_lim = 2, n_dp_aero = 3, n_mp_aero = 4, n_c_aero = 5, n_pres = 6    ! Indices for each parameter
    logical, save :: set_params(n_params) = .false.                     ! Logical telling if a parameter has been stored at least once
    logical, save :: printed_info(n_params) = .false.                   ! Logical telling if some info on a parameter has been printed
    character(len=200), dimension(n_params) :: str_def
    character(len=200) :: str_err_aero
    
    ! Info / error messages
    str_def(n_mp_aero)(:) = 'No input particle mass stored: Calculating masses from diameters assuming default density'
    str_def(n_pres)(:) = 'No input pressure stored: Using the default value'
    str_err_aero = 'Cannot access aerosol parameters: All parameters have not been stored'
    
    ! Input data to store
    
    if (present(naero_in)) then
        
        if (naero_in .ne. naero_saved) then
            
            naero_saved = naero_in
            
            if (allocated(dp_aero_lim_saved)) deallocate(dp_aero_lim_saved)
            if (allocated(dp_aero_saved)) deallocate(dp_aero_saved)
            if (allocated(mp_aero_saved)) deallocate(mp_aero_saved)
            if (allocated(c_aero_saved)) deallocate(c_aero_saved)
            
            allocate(dp_aero_lim_saved(naero_saved+1))
            allocate(dp_aero_saved(naero_saved))
            allocate(mp_aero_saved(naero_saved))
            allocate(c_aero_saved(naero_saved))
            
        end if
        
        set_params(n_naero) = .true.
        
    end if
    
    if (present(dp_aero_lim_in)) then
        dp_aero_lim_saved = dp_aero_lim_in
        set_params(n_dp_aero_lim) = .true.
    end if
    
    if (present(dp_aero_in)) then
        dp_aero_saved = dp_aero_in
        set_params(n_dp_aero) = .true.
    end if
    
    if (present(mp_aero_in)) then
        mp_aero_saved = mp_aero_in
        set_params(n_mp_aero) = .true.
    end if
    
    if (present(c_aero_in)) then
        c_aero_saved = c_aero_in
        set_params(n_c_aero) = .true.
    end if
    
    if (present(pres_in)) then
        pres_saved = pres_in
        set_params(n_pres) = .true.
    end if
    
    ! Output data to return
    
    if (present(naero)) then
        if (.not. set_params(n_naero)) then
            write(*,*) str_err_aero
            stop
        end if
        naero = naero_saved
    end if
    
    if (present(dp_aero_lim)) then
        if (.not. set_params(n_dp_aero_lim)) then
            write(*,*) str_err_aero
            stop
        end if
        dp_aero_lim = dp_aero_lim_saved
    end if
    
    if (present(dp_aero)) then
        if (.not. set_params(n_dp_aero)) then
            write(*,*) str_err_aero
            stop
        end if
        dp_aero = dp_aero_saved
    end if
    
    if (present(mp_aero)) then
        if (.not. set_params(n_mp_aero)) then
            if (.not. set_params(n_dp_aero)) then
                write(*,*) str_err_aero
                stop
            else
                ! Assess the mass in case it is not given as input
                mp_aero_saved = rhop_def*pi/6.d0*dp_aero_saved**3.d0
                if (.not. printed_info(n_mp_aero)) then
                    write(*,*) str_def(n_mp_aero)
                    printed_info(n_mp_aero) = .true.
                end if
            end if
        end if
        mp_aero = mp_aero_saved
    end if
    
    if (present(c_aero)) then
        if (.not. set_params(n_c_aero)) then
            write(*,*) str_err_aero
            stop
        end if
        c_aero = c_aero_saved
    end if
    
    if (present(pres)) then
        if (.not. printed_info(n_pres)) then
            if (.not. set_params(n_pres)) then
                write(*,*) str_def(n_pres)
            end if
            ! No messages will be printed after the first check
            printed_info(n_pres) = .true.
        end if
        pres = pres_saved
    end if
    
end subroutine get_ambient_params

!--------------------------------------------------------------------------------------------------
! Test for including a user-defined cluster sink routine
!--------------------------------------------------------------------------------------------------

! subroutine input_sink_test(insink)
    
    ! implicit none
    
    ! ! Output
    ! real(kind(1.d0)), intent(out) :: insink(:)                                  ! Additional sink for each cluster (1/s)
    
    ! !insink = 0.d0
    ! insink = 1.d-2
    
! end subroutine input_sink_test

!--------------------------------------------------------------------------------------------------
! Limits of the size bins for size-classifying the output cluster concentrations
!--------------------------------------------------------------------------------------------------

subroutine get_bin_limits(bin_limits)

    implicit none
    
    real(kind(1.d0)) :: bin_limits(nbins+1)     ! Vector of the bin limits (m)

    ! Bin limits according to mobility diameter
    bin_limits = (/1.05d0, 1.28d0, 1.73d0, 2.59d0, 4.27d0, 6.36d0/)*1.d-9 ! m

end subroutine get_bin_limits

!--------------------------------------------------------------------------------------------------
! Parameters related to size-classification, including the indices of clusters for each size bin,
! numbers of clusters in bins, ...
!--------------------------------------------------------------------------------------------------

subroutine group_size_bins(grouped_indices,nclust_per_bin)

    implicit none
    
    integer :: grouped_indices(0:nbins,nclust_max)          ! Indices of the clusters belonging to each bin
    integer :: nclust_per_bin(0:nbins)                      ! Total number of clusters belonging to each bin
    real(kind(1.d0)) :: mi, ri, ivalue
    real(kind(1.d0)) :: bin_limits(nbins+1)                 ! Vector of the bin limits (m)
    integer :: i, j, nclust_syst, indices(nclust_max,n_mol_types)
    logical :: lfound

    grouped_indices = 0
    nclust_per_bin =  0
    lfound = .false.
    
    if (.not. small_set_mode) then              ! The conditional is needed to compile also in the small set mode

        if(n_mol_types .gt. 1) then
            call get_molecule_numbers(nclust_syst,indices)
        else
            nclust_syst = nclust_max
            indices = reshape((/((i,i=1,nclust_syst),j=1,n_mol_types)/), shape(indices))
        end if

        call get_bin_limits(bin_limits)
        do i=1,nclust_syst
            if (sum(indices(i,:)) .eq. 1) then
                cycle                           ! Don't add monomers to the size bins
            end if
            call get_masses_and_radii(mi,ri,indices(i,:))
            ivalue = 2.d0*ri + 0.3d0*1.d-9      ! Assuming d_mob = d_mass + 0.3 nm
            if (ivalue .ge. bin_limits(1)) then
                lfound = .false.
                do j=1,nbins
                    if ((ivalue .ge. bin_limits(j)) .and. (ivalue .lt. bin_limits(j+1))) then
                        nclust_per_bin(j) = nclust_per_bin(j) + 1
                        grouped_indices(j,nclust_per_bin(j)) = i
                        lfound = .true.
                        exit
                    end if
                end do
                if (.not. lfound) then
                    write(*,*) 'Could not group cluster ', i, ', molecules: ', indices(i,:)
                    stop
                end if
            else
                nclust_per_bin(0) = nclust_per_bin(0) + 1
                grouped_indices(0,nclust_per_bin(0)) = i
                !write(*,*) 'Cluster below the lowest bin limit: ' , i, ', molecules: ', indices(i,:)
            end if
        end do
    
    end if

end subroutine group_size_bins


end module acdc_simulation_setup
