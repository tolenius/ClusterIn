module acdc_aerosol_parameters

implicit none

! Parameters and subroutines related to aerosol distributions and properties

real(kind(1.d0)), parameter :: ctot_aero = 1.d9             ! Total aerosol concentration in "average" atmosphere (m^-3)

integer, parameter :: naero_fixed = 10                      ! Number of aerosol bins
real(kind(1.d0)), parameter :: dp_min=1.d0*1.d-9            ! Aerosol size range (m)
real(kind(1.d0)), parameter :: dp_max=10.d0*1.d-6

real(kind(1.d0)), parameter :: dmed_aero = 100.d0*1.d-9     ! Median diameter (m)
real(kind(1.d0)), parameter :: sigma_aero = 2.d0            ! Geometric standard deviation
!real(kind(1.d0)), parameter :: rho_aero = 1500.d0           ! Particle density (kg m^-3)
!real(kind(1.d0)), parameter :: cs_ref = 3.95d-3             ! cs_ref corresponding to 100 nm, sigma = 2

real(kind(1.d0)), parameter :: pi=4.d0*atan(1.d0)


contains

!--------------------------------------------------------------------------------------------------
! Scavenging sink for each cluster onto each aerosol size bin,
! based on the exact aerosol size distribution and given coagulation coefficient formulae
!--------------------------------------------------------------------------------------------------

subroutine get_cs_aero(nclust,naero,dp_aero,mp_aero,c_aero,dp_clust,mp_clust,temperature,pressure,cs_per_clust_aero)

    implicit none
    
    integer, intent(in) :: nclust, naero                                        ! Numbers of clusters and aerosol bins
    real(kind(1.d0)), intent(in) :: dp_aero(naero), mp_aero(naero), c_aero(naero)   ! Aerosol diameters (m), masses (kg) and concentrations (m^-3)
    real(kind(1.d0)), intent(in) :: dp_clust(nclust), mp_clust(nclust)          ! Cluster diameters (m) and masses (kg)
    real(kind(1.d0)), intent(in) :: temperature                                 ! Temperature (K)
    real(kind(1.d0)), intent(in) :: pressure                                    ! Pressure (Pa)
    real(kind(1.d0)), intent(out) :: cs_per_clust_aero(naero,nclust)            ! Sink for each cluster onto each aerosol bin (1/s)
    real(kind(1.d0)) :: dp_all(nclust+naero), m_all(nclust+naero)
    real(kind(1.d0)) :: Cc(nclust+naero), Diff(nclust+naero), veloc(nclust+naero), mfp(nclust+naero), g(nclust+naero)
    real(kind(1.d0)) :: mu, lambda, fs_corr
    integer :: i, j
    real(kind(1.d0)), parameter :: m_air=28.97d0*1.d-3, Rg=8.3145d0, kB=1.3806504d-23
    
    ! Combine the cluster and aerosol diameters and masses into single arrays
    dp_all(1:nclust)=dp_clust
    dp_all(nclust+1:nclust+naero)=dp_aero
    m_all(1:nclust)=mp_clust
    m_all(nclust+1:nclust+naero)=mp_aero
    
    ! Parameters for the Fuchs-Sutugin coagulation coefficients
    mu=1.8d-5*(temperature/298.d0)**0.85d0                                      ! Viscosity of air (kg/m/s) (ADCHAM)
    lambda=2.d0*mu/(pressure*(8.d0*m_air/(pi*Rg*temperature))**0.5d0)           ! Gas mean free path in air (m) (Seinfeld and Pandis, Eq. 9.6)
    
    Cc=1.d0+2.d0*lambda/dp_all*(1.257d0+0.4d0*exp(-1.1d0*dp_all/(2.d0*lambda))) ! Slip correction factor (Seinfeld and Pandis, Eq. 9.34)
    Diff=Cc*kB*temperature/(3.d0*pi*mu*dp_all)                                  ! Particle diffusivity (m^2/s)
    
    veloc=(8.d0*kB*temperature/(pi*m_all))**0.5d0                               ! Thermal velocity (m/s)
    mfp=8.d0*Diff/(pi*veloc)                                                    ! Mean free path (m) (Seinfeld and Pandis, Table 13.1)

    g=1.d0/(3.d0*dp_all*mfp)*((dp_all+mfp)**3.d0-(dp_all**2.d0+mfp**2.d0)**(3.d0/2.d0))-dp_all
    
    ! Calculate the coagulation sink of each cluster onto each aerosol bin
    ! NOTE that this includes the aerosol bin concentration, i.e. the unit is 1/s, not m^3/s
    cs_per_clust_aero = 0.d0
    
    do i = 1, nclust
        do j = 1, naero
        
            ! Seinfeld and Pandis, Table 13.1
            fs_corr=((dp_all(i)+dp_all(nclust+j))/(dp_all(i)+dp_all(nclust+j)+2.d0*(g(i)**2.d0+g(nclust+j)**2.d0)**0.5d0)&
                &+8.d0*(Diff(i)+Diff(nclust+j))/((veloc(i)**2.d0+veloc(nclust+j)**2.d0)**0.5d0&
                &*(dp_all(i)+dp_all(nclust+j))))**(-1.d0)
            
            cs_per_clust_aero(j,i) = 2.d0*pi*(dp_all(i)+dp_all(nclust+j))*(Diff(i)+Diff(nclust+j))*fs_corr*c_aero(j)
            
        end do
    end do

end subroutine get_cs_aero

!--------------------------------------------------------------------------------------------------
! Get a lognormal aerosol distribution for given parameters
!--------------------------------------------------------------------------------------------------

subroutine lognormal_distr(c_tot,dp_med,gsd,nbins_distr,c_distr,dp_distr)

    implicit none
    
    real(kind(1.d0)), intent(in) :: c_tot                                       ! Total aerosol concentration (m^-3)
    real(kind(1.d0)), intent(in) :: dp_med                                      ! Median diameter (m)
    real(kind(1.d0)), intent(in) :: gsd                                         ! Standard deviation
    integer, intent(in) :: nbins_distr                                          ! Number of aerosol size bins
    real(kind(1.d0)), intent(out) :: c_distr(nbins_distr)                       ! Size-binned aerosol concentration (m^-3)
    real(kind(1.d0)), intent(out) :: dp_distr(nbins_distr)                      ! Bin diameter (m)
    
    real(kind(1.d0)) :: dp_lim(nbins_distr+1)                                   ! Bin limits (m)
    real(kind(1.d0)) :: dlogdp(nbins_distr)                                     ! Bin width (m)
    
    ! Example values:
    !c_tot = 1.d3*1.d6
    !dp_med = 100.d0*1.d-9
    !gsd = 2.d0
    
    call dp_lim_aero(nbins_distr,dp_lim,dp_distr)
    
    dlogdp=log10(dp_lim(2:nbins_distr+1)/dp_lim(1:nbins_distr))
    
    ! Log-normal number concentration density assuming a single mode (m^-3)
    ! dc/dlogdp
    c_distr=c_tot/((2.d0*pi)**0.5d0*log10(gsd))*exp(-(log10(dp_distr)-log10(dp_med))**2.d0/(2.d0*(log10(gsd))**2.d0))
    ! c in bins
    c_distr=c_distr*dlogdp
    
end subroutine lognormal_distr

!--------------------------------------------------------------------------------------------------
! Get aerosol bin diameter limits and mean diameters
!--------------------------------------------------------------------------------------------------

subroutine dp_lim_aero(nbins_distr,dp_lim,dp_mean)

    implicit none
    
    integer, intent(in) :: nbins_distr                                          ! Number of aerosol size bins
    real(kind(1.d0)), intent(out) :: dp_lim(nbins_distr+1)                      ! Bin limits (m)
    real(kind(1.d0)), intent(out) :: dp_mean(nbins_distr)                       ! Bin geometric mean diameter (m)
    
    real(kind(1.d0)) :: real_tmp
    integer :: i
    
    dp_lim=0.d0
    real_tmp=log10(dp_max/dp_min)/nbins_distr
    dp_lim=(/ (dp_min*10.d0**(i*real_tmp), i=0,nbins_distr) /)
    
    ! Geometric mean diameter in each bin (m)
    dp_mean=(dp_lim(1:nbins_distr)*dp_lim(2:nbins_distr+1))**0.5d0
    
end subroutine dp_lim_aero

end module acdc_aerosol_parameters