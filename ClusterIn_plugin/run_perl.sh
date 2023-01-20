#!/bin/bash


###############################################################################
####################             Input section             ####################
###############################################################################

#### Set the chemical system and optionally some conditions ####

# Vapor names (as in the input files for Perl)
vapors=("A" "N")

# Include charged clusters
l_incl_ions=1

# Fixed T and RH values - comment out when using varying input T
# (RH works only together with fixed T)
#temperature=280
#rh=20


#### Options for dynamics and coupling of cluster dynamics to an aerosol dynamics model ####

# Full dynamic coupling including aerosol evaporation, cluster scavenging mass transfer and size-classified J
l_full_dynamic=1
# Faster compilation: Vapor concentrations always forced constant - use carefully!
# Constant concentrations can also be set within ACDC when l_const_vapor=0, but compilation may be slow for complex systems
l_const_vapor=0


#### Additional options for generating the cluster equations by the Perl code - comment out if not used ####

# Tag that will be available in the acdc_system file
perl_opt_add+=" --tag test"


#### Suffices of the cluster system files which will be generated ####

# Suffix to add to the file names
file_suffix="_test"
# Include also the vapor names in the suffix
l_add_vapor_suffix=1

# Suffix to add to the subroutine names within the files
#routine_suffix="_test"


#### Paths and names of the Perl input files ####

data_dir="./Perl_input"

g_file="HS_YesSymYesQuasiHCorr_Besel.txt"
dip_file="dipoles_YesQuasiHCorr_Besel.txt"

# By default, the cluster set file is assumed to be "input_XYZ_neutral[_neg_pos][_inp_suffix].inp", where
# X, Y and Z are vapor names and the inclusion of "_neg_pos" depends on if ions are included;
# if inp_suffix is not used, comment out
inp_suffix="_Besel"

perl_file="acdc_2022_11_24.pl"


# End of the input section
###############################################################################




######## Determine the Perl input and create the cluster system files ########


[ $l_full_dynamic -eq 1 ] && l_const_vapor=0
[ $l_const_vapor -eq 1 ] && echo "Warning: Vapor concentrations always forced constant - this cannot be changed anymore in ACDC"

# Create the Perl option string

perl_opt=" --fortran --e $data_dir/$g_file"

# Cluster set name (assumed to correspond to vapor names) and conditions

cluster_file="input_"
vapor_suffix="_"

if [ -n "$temperature" ]; then
    perl_opt+=" --temperature $temperature"
else
    perl_opt+=" --variable_temp"
fi

if [ -n "$rh" ]; then
    perl_opt+=" --rh $rh"
fi

for vapor in "${vapors[@]}"; do

    perl_opt+=" --cs_only 1$vapor,0"
    [ $l_const_vapor -eq 1 ] && perl_opt+=" --no_eq 1$vapor"
    
    cluster_file+="$vapor"
    vapor_suffix+="$vapor"
    
done

cluster_file+="_neutral"

if [ $l_incl_ions -eq 1 ]; then
    perl_opt+=" --variable_ion_source --dip $data_dir/$dip_file"
    cluster_file+="_neg_pos"
    vapor_suffix+="_ions"
else
    vapor_suffix+="_noions"
fi

[ -n "$inp_suffix" ] && cluster_file+="$inp_suffix"
cluster_file+=".inp"

perl_opt+=" --i $data_dir/$cluster_file"

###############################################################################

# Possible dynamic coupling to an aerosol dynamics model

if [ $l_full_dynamic -eq 1 ]; then
    echo "----------------------------------------------------------------------------"
    echo "Applying a full dynamic coupling to an aerosol dynamics model"
    echo "Warning: Ensure that the steady state assumption is TURNED OFF within ACDC"
    echo "----------------------------------------------------------------------------"
    
    perl_opt+=" --cs external --save_coag_per_clust --save_outgoing_clust"
else
    # Approximate the size-dependent cluster scavenging sink by an exponent form based on the H2SO4 sink
    perl_opt+=" --variable_cs --cs exp_loss --exp_loss_exponent -1.6 --exp_loss_ref_size 0.55 --save_outgoing"
fi

###############################################################################

# Possible additional options

[ -n "$perl_opt_add" ] && perl_opt+=" $perl_opt_add"

# Examples of alternative / additional options
# --scale_evap_factor -1.0
# --wl CLOUD4_simple

###############################################################################

# Add possible suffices to the output files

suffix_tmp=""
[ $l_add_vapor_suffix -eq 1 ] && suffix_tmp+="$vapor_suffix"
[ -n "$file_suffix" ] && suffix_tmp+="$file_suffix"
[ "$suffix_tmp" != "" ] && perl_opt+=" --append $suffix_tmp"

# Add potential suffices to the subroutine names
[ -n "$routine_suffix" ] && perl_opt+=" --append_to_subroutines $routine_suffix"

###############################################################################

# Generate the equations

printf "\nGenerating the cluster equations with Perl options\n${perl_opt//--/\n--}\n"

perl $( echo "$perl_file $perl_opt" )

###############################################################################

# Print also a file containing info on the vapor species
printf '%s\n' "${#vapors[@]}" "${vapors[@]}" > "cluster_chem_spec${file_suffix}.txt"
