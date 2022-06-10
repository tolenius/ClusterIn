FC = gfortran
FCFLAGS = -O3 -fcheck=bounds -finit-local-zero
#FCFLAGS = -O3 -fcheck=bounds -finit-local-zero -mcmodel=medium -pg


# Names of some cluster files (in case one wants to test different versions with different names)
run_sim = run_clusterin_example.f90
clust_sim = clusterin.f90
driver = driver_acdc_J.f90
eqs = acdc_equations_AN.f90
sys = acdc_system_AN.f90


# Executable
run: run_clusterin.o clusterin.o driver.o acdc_equations.o acdc_system.o acdc_simulation_setup.o solution_settings.o dvode.o \
     acdc_aerosol_parameters.o
	$(FC) $(FCFLAGS) $^ -o $@


# Main program
run_clusterin.o: $(run_sim) clusterin.o acdc_simulation_setup.o acdc_aerosol_parameters.o
	$(FC) $(FCFLAGS) -c $< -o $@

# Files specific for a given cluster system
clusterin.o: $(clust_sim) driver.o acdc_system.o acdc_simulation_setup.o
	$(FC) $(FCFLAGS) -c $< -o $@

driver.o: $(driver) acdc_system.o acdc_simulation_setup.o solution_settings.o
	$(FC) $(FCFLAGS) -c $< -o $@

acdc_equations.o: $(eqs) acdc_simulation_setup.o
	$(FC) $(FCFLAGS) -c $< -o $@

acdc_system.o: $(sys)
	$(FC) $(FCFLAGS) -c $< -o $@

acdc_simulation_setup.o: acdc_simulation_setup.f90 acdc_system.o acdc_aerosol_parameters.o
	$(FC) $(FCFLAGS) -c $< -o $@

# Common files for all cluster systems
acdc_aerosol_parameters.o: acdc_aerosol_parameters.f90
	$(FC) $(FCFLAGS) -c $< -o $@

solution_settings.o: solvers/solution_settings.f90
	$(FC) $(FCFLAGS) -c $< -o $@
    
dvode.o: solvers/dvode.f
	$(FC) $(FCFLAGS) -std=legacy -c $< -o $@


.PHONY: clean

clean:
	rm -f *.o *.mod run
