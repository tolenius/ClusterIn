module acdc_system

implicit none

character(len=*), parameter :: tag = 'test'			! tag for the simulation system

logical, parameter :: small_set_mode = .true.			! small set of clusters (default mode)

integer, parameter :: nclust = 65						! number of clusters, molecules and ions
integer, parameter :: neq = 245							! number of equations
integer, parameter :: nclust_max = 65, neq_max = 245	! same as nclust and neq in the small set mode

integer, parameter :: n_variables = 3
logical, parameter :: variable_temp = .true.
logical, parameter :: variable_cs = .false.
logical, parameter :: variable_ipr = .true.
real(kind(1.d0)), parameter :: rh = 0.d0					! RH in %
real(kind(1.d0)), parameter :: pw = 0.d0					! partial pressure of water in Pa

integer, parameter :: n_monomer_types = 4
integer, parameter :: n_charges = 3			! number of charging states
integer, parameter :: n1A = 1, n1N = 4, n1B = 25, n1P1N = 47, nneg = 64, npos = 65			! cluster indices for monomers and ions

integer, parameter :: nout_all(3) = (/72, 73, 74/), nout_neu = 72, nout_neg = 73, nout_pos = 74			! indices for outgoing fluxes

integer, parameter :: nclust_coag_saved = 65			! number of saved scavenged clusters
integer, parameter :: nclust_out_saved = 105				! number of saved outgrown clusters
integer, parameter :: nrange_coag_saved(2) = (/76, 140/)	! index range of saved scavenged clusters
integer, parameter :: nrange_out_saved(2) = (/141, 245/)	! index range of saved outgrown clusters

integer, parameter :: n_mol_types = 4
integer, parameter :: nmolA = 1, nmolB = 2, nmolP = 3, nmolN = 4			! molecule indices for the used species

integer, parameter :: n_1A_clusters = 3				! number molecules and clusters containing 1 A molecule

integer, parameter :: n_neutrals = 24			! number of neutral molecules and clusters
integer, parameter :: n_negatives = 23			! negative
integer, parameter :: n_positives = 18			! positive
integer, parameter :: nclust_nogen = 63			! number of clusters and molecules excluding generic ions

integer, parameter :: n_neutral_monomers = 2, neutral_monomers(2) = (/1, 4/)			! number and indices of neutral monomers
integer, parameter :: n_negative_monomers = 1, negative_monomers(1) = (/25/)			! negative
integer, parameter :: n_positive_monomers = 1, positive_monomers(1) = (/47/)			! positive

integer, parameter :: n_neutral_clusters = 22			! number of neutral clusters
integer, parameter :: n_negative_clusters = 21			! negative
integer, parameter :: n_positive_clusters = 16			! positive

real(kind(1.d0)), parameter :: mass_max = 691.72
real(kind(1.d0)), parameter :: diameter_max = 1.14
real(kind(1.d0)), parameter :: mob_diameter_max = 1.44
integer, parameter :: ij_ind_max(4) = (/6, 1, 1, 6/)		! maximum molecular content
integer, parameter :: n_bound = 11		! number of clusters at the system boundary
integer, parameter :: n_comp_out = 105		! number of clusters growing out in different compositions

integer, parameter :: nmols_out_neutral(1, 4) = reshape((/7, 0, 0, 6/),(/1, 4/))			! criteria for outgrowing neutrals
integer, parameter :: nmols_out_negative(1, 4) = reshape((/6, 1, 0, 3/),(/1, 4/))			! negatives
integer, parameter :: nmols_out_positive(1, 4) = reshape((/6, 0, 1, 7/),(/1, 4/))			! positives


contains

subroutine n_A_in_clusters(n_A)
	implicit none
	integer :: n_A(65)

	n_A = (/1, 2, 3, 0, 1, 2, 3, 1, 2, 3, &
		&4, 2, 3, 4, 5, 3, 4, 5, 6, 4, &
		&5, 6, 5, 6, 1, 2, 3, 4, 1, 2, &
		&3, 4, 2, 3, 4, 5, 3, 4, 5, 6, &
		&4, 5, 6, 5, 6, 6, 0, 1, 0, 1, &
		&2, 1, 2, 3, 2, 3, 4, 3, 4, 5, &
		&4, 5, 6, 0, 0/)

end subroutine n_A_in_clusters

subroutine clusters_with_1_A(cluster_numbers)
	implicit none
	integer :: cluster_numbers(3)

	cluster_numbers = (/1, 5, 8/)

end subroutine clusters_with_1_A

subroutine arrays(neutrals, negatives, positives)
	implicit none
	integer :: neutrals(24), negatives(23), positives(18)

	neutrals = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
		&21, 22, 23, 24/)
	negatives = (/25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, &
		&40, 41, 42, 43, 44, 45, 46, 64/)
	positives = (/47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, &
		&62, 63, 65/)

end subroutine arrays

subroutine cluster_arrays(neutral_clusters, negative_clusters, positive_clusters)
	implicit none
	integer :: neutral_clusters(22), negative_clusters(21), positive_clusters(16)

	neutral_clusters = (/2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
		&21, 22, 23, 24/)
	negative_clusters = (/26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, &
		&40, 41, 42, 43, 44, 45, 46/)
	positive_clusters = (/48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, &
		&62, 63/)

end subroutine cluster_arrays

subroutine get_charging_state(charging_state)
	implicit none
	integer :: charging_state(65)

	charging_state = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 2, 2, 2, 2, 2, 2, &
		&2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
		&2, 2, 2, 2, 2, 2, 3, 3, 3, 3, &
		&3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
		&3, 3, 3, 2, 3/)

end subroutine get_charging_state

subroutine get_mass(mass)
	implicit none
	real(kind(1.d0)) :: mass(65)

	mass = (/98.08, 196.16, 294.24, 17.04, 115.12, 213.20, 311.28, 132.16, 230.24, 328.32, &
		&426.40, 247.28, 345.36, 443.44, 541.52, 362.40, 460.48, 558.56, 656.64, 477.52, &
		&575.60, 673.68, 592.64, 690.72, 97.08, 195.16, 293.24, 391.32, 114.12, 212.20, &
		&310.28, 408.36, 229.24, 327.32, 425.40, 523.48, 344.36, 442.44, 540.52, 638.60, &
		&459.48, 557.56, 655.64, 574.60, 672.68, 689.72, 18.04, 116.12, 35.08, 133.16, &
		&231.24, 150.20, 248.28, 346.36, 265.32, 363.40, 461.48, 380.44, 478.52, 576.60, &
		&495.56, 593.64, 691.72, 32.00, 19.02/)

end subroutine get_mass

subroutine get_diameter(diameter)
	implicit none

	real(kind(1.d0)) :: diameter(65)

	 diameter = (/0.55, 0.70, 0.80, 0.43, 0.63, 0.75, 0.84, 0.69, 0.79, 0.87, &
		&0.94, 0.83, 0.91, 0.97, 1.03, 0.94, 1.00, 1.05, 1.10, 1.02, &
		&1.07, 1.12, 1.10, 1.14, 0.55, 0.70, 0.80, 0.88, 0.63, 0.75, &
		&0.84, 0.91, 0.79, 0.87, 0.94, 1.00, 0.90, 0.97, 1.03, 1.08, &
		&1.00, 1.05, 1.10, 1.07, 1.12, 1.14, 0.43, 0.63, 0.54, 0.69, &
		&0.79, 0.74, 0.83, 0.91, 0.87, 0.94, 1.00, 0.96, 1.02, 1.07, &
		&1.05, 1.10, 1.14, 0.45, 0.39/)	! dry value

end subroutine get_diameter

subroutine get_mob_diameter(mob_diameter)
	implicit none

	real(kind(1.d0)) :: mob_diameter(65)

	 mob_diameter = (/0.85, 1.00, 1.10, 0.73, 0.93, 1.05, 1.14, 0.99, 1.09, 1.17, &
		&1.24, 1.13, 1.21, 1.27, 1.33, 1.24, 1.30, 1.35, 1.40, 1.32, &
		&1.37, 1.42, 1.40, 1.44, 0.85, 1.00, 1.10, 1.18, 0.93, 1.05, &
		&1.14, 1.21, 1.09, 1.17, 1.24, 1.30, 1.20, 1.27, 1.33, 1.38, &
		&1.30, 1.35, 1.40, 1.37, 1.42, 1.44, 0.73, 0.93, 0.84, 0.99, &
		&1.09, 1.04, 1.13, 1.21, 1.17, 1.24, 1.30, 1.26, 1.32, 1.37, &
		&1.35, 1.40, 1.44, 0.75, 0.69/)	! dry value

end subroutine get_mob_diameter

subroutine cluster_names(clust)
	implicit none
	character(len=11), dimension(75) :: clust

	clust(1)(:) = '1A'
	clust(2)(:) = '2A'
	clust(3)(:) = '3A'
	clust(4)(:) = '1N'
	clust(5)(:) = '1A1N'
	clust(6)(:) = '2A1N'
	clust(7)(:) = '3A1N'
	clust(8)(:) = '1A2N'
	clust(9)(:) = '2A2N'
	clust(10)(:) = '3A2N'
	clust(11)(:) = '4A2N'
	clust(12)(:) = '2A3N'
	clust(13)(:) = '3A3N'
	clust(14)(:) = '4A3N'
	clust(15)(:) = '5A3N'
	clust(16)(:) = '3A4N'
	clust(17)(:) = '4A4N'
	clust(18)(:) = '5A4N'
	clust(19)(:) = '6A4N'
	clust(20)(:) = '4A5N'
	clust(21)(:) = '5A5N'
	clust(22)(:) = '6A5N'
	clust(23)(:) = '5A6N'
	clust(24)(:) = '6A6N'
	clust(25)(:) = '1B'
	clust(26)(:) = '1A1B'
	clust(27)(:) = '2A1B'
	clust(28)(:) = '3A1B'
	clust(29)(:) = '1B1N'
	clust(30)(:) = '1A1B1N'
	clust(31)(:) = '2A1B1N'
	clust(32)(:) = '3A1B1N'
	clust(33)(:) = '1A1B2N'
	clust(34)(:) = '2A1B2N'
	clust(35)(:) = '3A1B2N'
	clust(36)(:) = '4A1B2N'
	clust(37)(:) = '2A1B3N'
	clust(38)(:) = '3A1B3N'
	clust(39)(:) = '4A1B3N'
	clust(40)(:) = '5A1B3N'
	clust(41)(:) = '3A1B4N'
	clust(42)(:) = '4A1B4N'
	clust(43)(:) = '5A1B4N'
	clust(44)(:) = '4A1B5N'
	clust(45)(:) = '5A1B5N'
	clust(46)(:) = '5A1B6N'
	clust(47)(:) = '1P1N'
	clust(48)(:) = '1A1P1N'
	clust(49)(:) = '1P2N'
	clust(50)(:) = '1A1P2N'
	clust(51)(:) = '2A1P2N'
	clust(52)(:) = '1A1P3N'
	clust(53)(:) = '2A1P3N'
	clust(54)(:) = '3A1P3N'
	clust(55)(:) = '2A1P4N'
	clust(56)(:) = '3A1P4N'
	clust(57)(:) = '4A1P4N'
	clust(58)(:) = '3A1P5N'
	clust(59)(:) = '4A1P5N'
	clust(60)(:) = '5A1P5N'
	clust(61)(:) = '4A1P6N'
	clust(62)(:) = '5A1P6N'
	clust(63)(:) = '6A1P6N'
	clust(64)(:) = 'neg'
	clust(65)(:) = 'pos'
	clust(66)(:) = 'source'
	clust(67)(:) = 'coag'
	clust(68)(:) = 'wall'
	clust(69)(:) = 'dil'
	clust(70)(:) = 'insink'
	clust(71)(:) = 'rec'
	clust(72)(:) = 'out_neu'
	clust(73)(:) = 'out_neg'
	clust(74)(:) = 'out_pos'
	clust(75)(:) = 'bound'

	! Additional flux elements not listed here:
	! Elements from 76 to 140: coagulation loss of each cluster
	! Elements from 141 to 245: each outgrown cluster composition

end subroutine cluster_names

subroutine monomer_names(clust_mon)
	implicit none
	character(len=11), dimension(4) :: clust_mon

	clust_mon(1)(:) = '1A'
	clust_mon(2)(:) = '1N'
	clust_mon(3)(:) = '1B'
	clust_mon(4)(:) = '1P1N'

end subroutine monomer_names

subroutine molecule_names(labels)
	implicit none
	character(len=11), dimension(4) :: labels

	labels(1)(:) = 'A'
	labels(2)(:) = 'B'
	labels(3)(:) = 'P'
	labels(4)(:) = 'N'

end subroutine molecule_names

subroutine monomer_indices(n_monomers)
	implicit none
	integer :: n_monomers(4)

	n_monomers = (/1, 4, 25, 47/)

end subroutine monomer_indices

subroutine get_bound(bound_clusters,nmols_bound)
	implicit none
	integer :: bound_clusters(11), nmols_bound(11,4)

	nmols_bound(1,:) = (/6, 0, 0, 4/)
	nmols_bound(2,:) = (/6, 0, 0, 5/)
	nmols_bound(3,:) = (/5, 0, 0, 6/)
	nmols_bound(4,:) = (/6, 0, 0, 6/)
	nmols_bound(5,:) = (/5, 1, 0, 3/)
	nmols_bound(6,:) = (/5, 1, 0, 4/)
	nmols_bound(7,:) = (/5, 1, 0, 5/)
	nmols_bound(8,:) = (/5, 1, 0, 6/)
	nmols_bound(9,:) = (/4, 0, 1, 6/)
	nmols_bound(10,:) = (/5, 0, 1, 6/)
	nmols_bound(11,:) = (/6, 0, 1, 6/)

	bound_clusters = (/19, 22, 23, 24, 40, 43, 45, 46, 61, 62, 63/)

end subroutine get_bound

subroutine get_diameter_out(diameter_out)
	implicit none
	real(kind(1.d0)) :: diameter_out(105)

	diameter_out = (/1.18, 1.12, 1.14, 1.16, 1.18, 1.22, 1.17, 1.19, 1.20, 1.22, &
		&1.26, 1.21, 1.23, 1.24, 1.26, 1.16, 1.20, 1.20, 1.20, 1.24, &
		&1.24, 1.24, 1.28, 1.27, 1.28, 1.22, 1.22, 1.18, 1.22, 1.26, &
		&1.26, 1.26, 1.29, 1.29, 1.29, 1.29, 1.31, 1.32, 1.28, 1.29, &
		&1.31, 1.32, 1.32, 1.24, 1.27, 1.27, 1.20, 1.24, 1.27, 1.31, &
		&1.31, 1.31, 1.34, 1.34, 1.34, 1.34, 1.36, 1.37, 1.33, 1.34, &
		&1.36, 1.37, 1.37, 1.29, 1.32, 1.32, 1.25, 1.29, 1.32, 1.35, &
		&1.35, 1.35, 1.38, 1.38, 1.38, 1.39, 1.40, 1.41, 1.26, 1.37, &
		&1.39, 1.40, 1.41, 1.41, 1.34, 1.37, 1.37, 1.30, 1.34, 1.37, &
		&1.40, 1.40, 1.40, 1.43, 1.42, 1.43, 1.38, 1.41, 1.41, 1.35, &
		&1.38, 1.41, 1.44, 1.44, 1.44/)	! dry value

end subroutine get_diameter_out

subroutine get_nmols_out(nmols_out)
	implicit none
	integer :: nmols_out(105,4)

	nmols_out(1,:) = (/7, 0, 0, 6/)
	nmols_out(2,:) = (/6, 1, 0, 3/)
	nmols_out(3,:) = (/6, 1, 0, 4/)
	nmols_out(4,:) = (/6, 1, 0, 5/)
	nmols_out(5,:) = (/6, 1, 0, 6/)
	nmols_out(6,:) = (/8, 0, 0, 6/)
	nmols_out(7,:) = (/7, 1, 0, 3/)
	nmols_out(8,:) = (/7, 1, 0, 4/)
	nmols_out(9,:) = (/7, 1, 0, 5/)
	nmols_out(10,:) = (/7, 1, 0, 6/)
	nmols_out(11,:) = (/9, 0, 0, 6/)
	nmols_out(12,:) = (/8, 1, 0, 3/)
	nmols_out(13,:) = (/8, 1, 0, 4/)
	nmols_out(14,:) = (/8, 1, 0, 5/)
	nmols_out(15,:) = (/8, 1, 0, 6/)
	nmols_out(16,:) = (/6, 0, 1, 7/)
	nmols_out(17,:) = (/7, 0, 0, 7/)
	nmols_out(18,:) = (/6, 1, 0, 7/)
	nmols_out(19,:) = (/7, 0, 1, 7/)
	nmols_out(20,:) = (/8, 0, 0, 7/)
	nmols_out(21,:) = (/7, 1, 0, 7/)
	nmols_out(22,:) = (/8, 0, 1, 7/)
	nmols_out(23,:) = (/9, 0, 0, 7/)
	nmols_out(24,:) = (/8, 1, 0, 7/)
	nmols_out(25,:) = (/9, 0, 1, 7/)
	nmols_out(26,:) = (/7, 0, 0, 8/)
	nmols_out(27,:) = (/6, 1, 0, 8/)
	nmols_out(28,:) = (/6, 0, 1, 8/)
	nmols_out(29,:) = (/7, 0, 1, 8/)
	nmols_out(30,:) = (/8, 0, 0, 8/)
	nmols_out(31,:) = (/7, 1, 0, 8/)
	nmols_out(32,:) = (/8, 0, 1, 8/)
	nmols_out(33,:) = (/9, 0, 0, 8/)
	nmols_out(34,:) = (/8, 1, 0, 8/)
	nmols_out(35,:) = (/9, 0, 1, 8/)
	nmols_out(36,:) = (/10, 0, 0, 6/)
	nmols_out(37,:) = (/10, 0, 0, 7/)
	nmols_out(38,:) = (/10, 0, 0, 8/)
	nmols_out(39,:) = (/9, 1, 0, 5/)
	nmols_out(40,:) = (/9, 1, 0, 6/)
	nmols_out(41,:) = (/9, 1, 0, 7/)
	nmols_out(42,:) = (/9, 1, 0, 8/)
	nmols_out(43,:) = (/10, 0, 1, 8/)
	nmols_out(44,:) = (/7, 0, 0, 9/)
	nmols_out(45,:) = (/8, 0, 0, 9/)
	nmols_out(46,:) = (/7, 1, 0, 9/)
	nmols_out(47,:) = (/6, 0, 1, 9/)
	nmols_out(48,:) = (/7, 0, 1, 9/)
	nmols_out(49,:) = (/8, 0, 1, 9/)
	nmols_out(50,:) = (/9, 0, 0, 9/)
	nmols_out(51,:) = (/8, 1, 0, 9/)
	nmols_out(52,:) = (/9, 0, 1, 9/)
	nmols_out(53,:) = (/10, 0, 0, 9/)
	nmols_out(54,:) = (/9, 1, 0, 9/)
	nmols_out(55,:) = (/10, 0, 1, 9/)
	nmols_out(56,:) = (/11, 0, 0, 7/)
	nmols_out(57,:) = (/11, 0, 0, 8/)
	nmols_out(58,:) = (/11, 0, 0, 9/)
	nmols_out(59,:) = (/10, 1, 0, 6/)
	nmols_out(60,:) = (/10, 1, 0, 7/)
	nmols_out(61,:) = (/10, 1, 0, 8/)
	nmols_out(62,:) = (/10, 1, 0, 9/)
	nmols_out(63,:) = (/11, 0, 1, 9/)
	nmols_out(64,:) = (/8, 0, 0, 10/)
	nmols_out(65,:) = (/9, 0, 0, 10/)
	nmols_out(66,:) = (/8, 1, 0, 10/)
	nmols_out(67,:) = (/7, 0, 1, 10/)
	nmols_out(68,:) = (/8, 0, 1, 10/)
	nmols_out(69,:) = (/9, 0, 1, 10/)
	nmols_out(70,:) = (/10, 0, 0, 10/)
	nmols_out(71,:) = (/9, 1, 0, 10/)
	nmols_out(72,:) = (/10, 0, 1, 10/)
	nmols_out(73,:) = (/11, 0, 0, 10/)
	nmols_out(74,:) = (/10, 1, 0, 10/)
	nmols_out(75,:) = (/11, 0, 1, 10/)
	nmols_out(76,:) = (/12, 0, 0, 8/)
	nmols_out(77,:) = (/12, 0, 0, 9/)
	nmols_out(78,:) = (/12, 0, 0, 10/)
	nmols_out(79,:) = (/9, 1, 0, 4/)
	nmols_out(80,:) = (/11, 1, 0, 7/)
	nmols_out(81,:) = (/11, 1, 0, 8/)
	nmols_out(82,:) = (/11, 1, 0, 9/)
	nmols_out(83,:) = (/11, 1, 0, 10/)
	nmols_out(84,:) = (/12, 0, 1, 10/)
	nmols_out(85,:) = (/9, 0, 0, 11/)
	nmols_out(86,:) = (/10, 0, 0, 11/)
	nmols_out(87,:) = (/9, 1, 0, 11/)
	nmols_out(88,:) = (/8, 0, 1, 11/)
	nmols_out(89,:) = (/9, 0, 1, 11/)
	nmols_out(90,:) = (/10, 0, 1, 11/)
	nmols_out(91,:) = (/11, 0, 0, 11/)
	nmols_out(92,:) = (/10, 1, 0, 11/)
	nmols_out(93,:) = (/11, 0, 1, 11/)
	nmols_out(94,:) = (/12, 0, 0, 11/)
	nmols_out(95,:) = (/11, 1, 0, 11/)
	nmols_out(96,:) = (/12, 0, 1, 11/)
	nmols_out(97,:) = (/10, 0, 0, 12/)
	nmols_out(98,:) = (/11, 0, 0, 12/)
	nmols_out(99,:) = (/10, 1, 0, 12/)
	nmols_out(100,:) = (/9, 0, 1, 12/)
	nmols_out(101,:) = (/10, 0, 1, 12/)
	nmols_out(102,:) = (/11, 0, 1, 12/)
	nmols_out(103,:) = (/12, 0, 0, 12/)
	nmols_out(104,:) = (/11, 1, 0, 12/)
	nmols_out(105,:) = (/12, 0, 1, 12/)

end subroutine get_nmols_out

subroutine molecule_names_nocharge(labels_nocharge)
	implicit none
	character(len=11), dimension(2) :: labels_nocharge

	labels_nocharge(1)(:) = 'A'
	labels_nocharge(2)(:) = 'N'

end subroutine molecule_names_nocharge

subroutine get_nmols_nocharge(nmols_nocharge_clust)
	implicit none
	integer :: nmols_nocharge_clust(65,2)

	nmols_nocharge_clust = 0

	nmols_nocharge_clust(1,:) = (/1, 0/)
	nmols_nocharge_clust(2,:) = (/2, 0/)
	nmols_nocharge_clust(3,:) = (/3, 0/)
	nmols_nocharge_clust(4,:) = (/0, 1/)
	nmols_nocharge_clust(5,:) = (/1, 1/)
	nmols_nocharge_clust(6,:) = (/2, 1/)
	nmols_nocharge_clust(7,:) = (/3, 1/)
	nmols_nocharge_clust(8,:) = (/1, 2/)
	nmols_nocharge_clust(9,:) = (/2, 2/)
	nmols_nocharge_clust(10,:) = (/3, 2/)
	nmols_nocharge_clust(11,:) = (/4, 2/)
	nmols_nocharge_clust(12,:) = (/2, 3/)
	nmols_nocharge_clust(13,:) = (/3, 3/)
	nmols_nocharge_clust(14,:) = (/4, 3/)
	nmols_nocharge_clust(15,:) = (/5, 3/)
	nmols_nocharge_clust(16,:) = (/3, 4/)
	nmols_nocharge_clust(17,:) = (/4, 4/)
	nmols_nocharge_clust(18,:) = (/5, 4/)
	nmols_nocharge_clust(19,:) = (/6, 4/)
	nmols_nocharge_clust(20,:) = (/4, 5/)
	nmols_nocharge_clust(21,:) = (/5, 5/)
	nmols_nocharge_clust(22,:) = (/6, 5/)
	nmols_nocharge_clust(23,:) = (/5, 6/)
	nmols_nocharge_clust(24,:) = (/6, 6/)
	nmols_nocharge_clust(25,:) = (/1, 0/)
	nmols_nocharge_clust(26,:) = (/2, 0/)
	nmols_nocharge_clust(27,:) = (/3, 0/)
	nmols_nocharge_clust(28,:) = (/4, 0/)
	nmols_nocharge_clust(29,:) = (/1, 1/)
	nmols_nocharge_clust(30,:) = (/2, 1/)
	nmols_nocharge_clust(31,:) = (/3, 1/)
	nmols_nocharge_clust(32,:) = (/4, 1/)
	nmols_nocharge_clust(33,:) = (/2, 2/)
	nmols_nocharge_clust(34,:) = (/3, 2/)
	nmols_nocharge_clust(35,:) = (/4, 2/)
	nmols_nocharge_clust(36,:) = (/5, 2/)
	nmols_nocharge_clust(37,:) = (/3, 3/)
	nmols_nocharge_clust(38,:) = (/4, 3/)
	nmols_nocharge_clust(39,:) = (/5, 3/)
	nmols_nocharge_clust(40,:) = (/6, 3/)
	nmols_nocharge_clust(41,:) = (/4, 4/)
	nmols_nocharge_clust(42,:) = (/5, 4/)
	nmols_nocharge_clust(43,:) = (/6, 4/)
	nmols_nocharge_clust(44,:) = (/5, 5/)
	nmols_nocharge_clust(45,:) = (/6, 5/)
	nmols_nocharge_clust(46,:) = (/6, 6/)
	nmols_nocharge_clust(47,:) = (/0, 1/)
	nmols_nocharge_clust(48,:) = (/1, 1/)
	nmols_nocharge_clust(49,:) = (/0, 2/)
	nmols_nocharge_clust(50,:) = (/1, 2/)
	nmols_nocharge_clust(51,:) = (/2, 2/)
	nmols_nocharge_clust(52,:) = (/1, 3/)
	nmols_nocharge_clust(53,:) = (/2, 3/)
	nmols_nocharge_clust(54,:) = (/3, 3/)
	nmols_nocharge_clust(55,:) = (/2, 4/)
	nmols_nocharge_clust(56,:) = (/3, 4/)
	nmols_nocharge_clust(57,:) = (/4, 4/)
	nmols_nocharge_clust(58,:) = (/3, 5/)
	nmols_nocharge_clust(59,:) = (/4, 5/)
	nmols_nocharge_clust(60,:) = (/5, 5/)
	nmols_nocharge_clust(61,:) = (/4, 6/)
	nmols_nocharge_clust(62,:) = (/5, 6/)
	nmols_nocharge_clust(63,:) = (/6, 6/)

end subroutine get_nmols_nocharge

subroutine get_nmols_nocharge_out(nmols_nocharge_out)
	implicit none
	integer :: nmols_nocharge_out(105,2)

	nmols_nocharge_out = 0

	nmols_nocharge_out(1,:) = (/7, 6/)
	nmols_nocharge_out(2,:) = (/7, 3/)
	nmols_nocharge_out(3,:) = (/7, 4/)
	nmols_nocharge_out(4,:) = (/7, 5/)
	nmols_nocharge_out(5,:) = (/7, 6/)
	nmols_nocharge_out(6,:) = (/8, 6/)
	nmols_nocharge_out(7,:) = (/8, 3/)
	nmols_nocharge_out(8,:) = (/8, 4/)
	nmols_nocharge_out(9,:) = (/8, 5/)
	nmols_nocharge_out(10,:) = (/8, 6/)
	nmols_nocharge_out(11,:) = (/9, 6/)
	nmols_nocharge_out(12,:) = (/9, 3/)
	nmols_nocharge_out(13,:) = (/9, 4/)
	nmols_nocharge_out(14,:) = (/9, 5/)
	nmols_nocharge_out(15,:) = (/9, 6/)
	nmols_nocharge_out(16,:) = (/6, 7/)
	nmols_nocharge_out(17,:) = (/7, 7/)
	nmols_nocharge_out(18,:) = (/7, 7/)
	nmols_nocharge_out(19,:) = (/7, 7/)
	nmols_nocharge_out(20,:) = (/8, 7/)
	nmols_nocharge_out(21,:) = (/8, 7/)
	nmols_nocharge_out(22,:) = (/8, 7/)
	nmols_nocharge_out(23,:) = (/9, 7/)
	nmols_nocharge_out(24,:) = (/9, 7/)
	nmols_nocharge_out(25,:) = (/9, 7/)
	nmols_nocharge_out(26,:) = (/7, 8/)
	nmols_nocharge_out(27,:) = (/7, 8/)
	nmols_nocharge_out(28,:) = (/6, 8/)
	nmols_nocharge_out(29,:) = (/7, 8/)
	nmols_nocharge_out(30,:) = (/8, 8/)
	nmols_nocharge_out(31,:) = (/8, 8/)
	nmols_nocharge_out(32,:) = (/8, 8/)
	nmols_nocharge_out(33,:) = (/9, 8/)
	nmols_nocharge_out(34,:) = (/9, 8/)
	nmols_nocharge_out(35,:) = (/9, 8/)
	nmols_nocharge_out(36,:) = (/10, 6/)
	nmols_nocharge_out(37,:) = (/10, 7/)
	nmols_nocharge_out(38,:) = (/10, 8/)
	nmols_nocharge_out(39,:) = (/10, 5/)
	nmols_nocharge_out(40,:) = (/10, 6/)
	nmols_nocharge_out(41,:) = (/10, 7/)
	nmols_nocharge_out(42,:) = (/10, 8/)
	nmols_nocharge_out(43,:) = (/10, 8/)
	nmols_nocharge_out(44,:) = (/7, 9/)
	nmols_nocharge_out(45,:) = (/8, 9/)
	nmols_nocharge_out(46,:) = (/8, 9/)
	nmols_nocharge_out(47,:) = (/6, 9/)
	nmols_nocharge_out(48,:) = (/7, 9/)
	nmols_nocharge_out(49,:) = (/8, 9/)
	nmols_nocharge_out(50,:) = (/9, 9/)
	nmols_nocharge_out(51,:) = (/9, 9/)
	nmols_nocharge_out(52,:) = (/9, 9/)
	nmols_nocharge_out(53,:) = (/10, 9/)
	nmols_nocharge_out(54,:) = (/10, 9/)
	nmols_nocharge_out(55,:) = (/10, 9/)
	nmols_nocharge_out(56,:) = (/11, 7/)
	nmols_nocharge_out(57,:) = (/11, 8/)
	nmols_nocharge_out(58,:) = (/11, 9/)
	nmols_nocharge_out(59,:) = (/11, 6/)
	nmols_nocharge_out(60,:) = (/11, 7/)
	nmols_nocharge_out(61,:) = (/11, 8/)
	nmols_nocharge_out(62,:) = (/11, 9/)
	nmols_nocharge_out(63,:) = (/11, 9/)
	nmols_nocharge_out(64,:) = (/8, 10/)
	nmols_nocharge_out(65,:) = (/9, 10/)
	nmols_nocharge_out(66,:) = (/9, 10/)
	nmols_nocharge_out(67,:) = (/7, 10/)
	nmols_nocharge_out(68,:) = (/8, 10/)
	nmols_nocharge_out(69,:) = (/9, 10/)
	nmols_nocharge_out(70,:) = (/10, 10/)
	nmols_nocharge_out(71,:) = (/10, 10/)
	nmols_nocharge_out(72,:) = (/10, 10/)
	nmols_nocharge_out(73,:) = (/11, 10/)
	nmols_nocharge_out(74,:) = (/11, 10/)
	nmols_nocharge_out(75,:) = (/11, 10/)
	nmols_nocharge_out(76,:) = (/12, 8/)
	nmols_nocharge_out(77,:) = (/12, 9/)
	nmols_nocharge_out(78,:) = (/12, 10/)
	nmols_nocharge_out(79,:) = (/10, 4/)
	nmols_nocharge_out(80,:) = (/12, 7/)
	nmols_nocharge_out(81,:) = (/12, 8/)
	nmols_nocharge_out(82,:) = (/12, 9/)
	nmols_nocharge_out(83,:) = (/12, 10/)
	nmols_nocharge_out(84,:) = (/12, 10/)
	nmols_nocharge_out(85,:) = (/9, 11/)
	nmols_nocharge_out(86,:) = (/10, 11/)
	nmols_nocharge_out(87,:) = (/10, 11/)
	nmols_nocharge_out(88,:) = (/8, 11/)
	nmols_nocharge_out(89,:) = (/9, 11/)
	nmols_nocharge_out(90,:) = (/10, 11/)
	nmols_nocharge_out(91,:) = (/11, 11/)
	nmols_nocharge_out(92,:) = (/11, 11/)
	nmols_nocharge_out(93,:) = (/11, 11/)
	nmols_nocharge_out(94,:) = (/12, 11/)
	nmols_nocharge_out(95,:) = (/12, 11/)
	nmols_nocharge_out(96,:) = (/12, 11/)
	nmols_nocharge_out(97,:) = (/10, 12/)
	nmols_nocharge_out(98,:) = (/11, 12/)
	nmols_nocharge_out(99,:) = (/11, 12/)
	nmols_nocharge_out(100,:) = (/9, 12/)
	nmols_nocharge_out(101,:) = (/10, 12/)
	nmols_nocharge_out(102,:) = (/11, 12/)
	nmols_nocharge_out(103,:) = (/12, 12/)
	nmols_nocharge_out(104,:) = (/12, 12/)
	nmols_nocharge_out(105,:) = (/12, 12/)

end subroutine get_nmols_nocharge_out


end module acdc_system

