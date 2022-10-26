#!/bin/bash


###############################################################################

# Default settings

l_verbose=0

# Directory where the files are located, if no input is given
dir="."
# Name of the subdirectories without the _1, _2, _3, ... suffix, if no input is given
subdir_dir="."
subdir_prefix="cluster_chem"
# Name of a Fortran include file that contains the "use" statements for the separate modules
inc_use_file="cluster_chem_use.inc"
# Name of the script that calls Perl to generate the cluster equation files
run_perl_file="run_perl.sh"


# Keyword which will be replaced by the system suffix (_1, _2, _3, ...)
suffix_key="chem_suffix"
# Delimiter used for splitting some array elements (see below)
delim=","


# All files that need to be copied
files_to_copy=( "$run_perl_file" "acdc_simulation_setup.f90" "clusterin.f90" "driver_acdc_J.f90" )


# All modules and subroutines that need to be renamed

# Search for the following patterns in the files to find the routine names
# Syntax: each element will be split at the delimiter, and the routine names are assumed to be the strings found in between the two parts
# Note that some characters must be escaped
strs_to_search=("module ${delim}" "subroutine ${delim}\(" "call ${delim}\(")

# Module / subroutine names that may not appear in the search but are included by default (i.e. are included in the Perl-generated files)
names_to_modify=("acdc_system" "feval" "jeval")
# Module / subroutine names that will not be modified (i.e. common routines that all sub-systems use)
names_to_exclude=("dvode" "get_cs_aero")

# Parameters to set in $run_perl_file which calls Perl
# Syntax: each element will be split at the delimiter, and the two parts give the parameter name and the parameter value
params_to_set=("file_suffix${delim}${suffix_key}" "routine_suffix${delim}${suffix_key}" "l_add_vapor_suffix${delim}0")


###############################################################################

# Get the command line input

help_str="\n
Generate separate subdirectories and files for different molecular cluster chemistries in order to include
several separate cluster dynamics systems simultaneously in an aerosol dynamics model

The cluster formation routines are copied into subdirectories sub_1, sub_2, sub_3, ...

Usage:
./copy_cluster_chem.sh -n number [-dir name] [-chemdir name] [-add name,name,...] [-exc name,name,...] [-v] [-h]

n          numbers of the separate chemical systems; e.g. 1-3 or 2,4
dir        name of the main directory that contains the files to be copied (optional, otherwise the current directory is used)
chemdir    name and/or path of the subdirectories (optional, otherwise a default name and the current directory are used)
add        additional files to copy, e.g. user-defined sink routines if they are not common for all chemistries (optional)
exc        routine / module names that are not renamed by _1, _2, ... in the additional files, i.e. those that are not chemistry-specific (optional)
v          verbose - it's good to use this option when you run the script for the first time to see that everything is ok
h          print these instructions
\n"

if [ $# -eq 0 ]; then
    printf "$help_str"
    exit 1
fi

while [ $# -gt 0 ]; do
    case $1 in
        -n) shift && nchem_str=$1;;
        -dir) shift && dir=$1;;
        -chemdir) shift && subdir_str=$1;;
        -add | -exc) str_tmp=$1
              shift
              readarray -t add_array_tmp < <( echo $1 | perl -pe 's/,\s*/\n/g' )
              if [ "$str_tmp" == "-add" ]; then
                  array_tmp=( "${files_to_copy[@]}" "${add_array_tmp[@]}" )
              elif [ "$str_tmp" == "-exc" ]; then
                  array_tmp=( "${names_to_exclude[@]}" "${add_array_tmp[@]}" )
              fi
              # Remove duplicates
              eval array_tmp=( $( printf "%q\n" "${array_tmp[@]}" | sort -u ) )
              if [ "$str_tmp" == "-add" ]; then
                  files_to_copy=( "${array_tmp[@]}" )
              elif [ "$str_tmp" == "-exc" ]; then
                  names_to_exclude=( "${array_tmp[@]}" )
              fi
              ;;
        -v) l_verbose=1;;
        -h | *) printf "$help_str"; exit 1;;
    esac
    shift
done

if [ -z "$nchem_str" ]; then
    printf "$help_str"
    exit 1
fi


###############################################################################

printf "\nCopying files for separate chemical system(s) $nchem_str\n\n"

readarray -t nchem_arr < <( echo "$nchem_str" | perl -pe 's/\s*//g' | perl -pe 's/-|,/\n/g' )
if [ "$nchem_str" != "${nchem_str//-/}" ]; then
    nchem_arr=( $( seq ${nchem_arr[0]} ${nchem_arr[1]} ) )
fi

# Get the names of all modules and subroutines in order to rename them for each chemical system

# Note: directory names are here without the last "/"
dir="${dir%/}"

# Add the possible directory path to the file names
files_to_copy=( "${files_to_copy[@]/#/$dir/}" )

# Add the possible directory path to the subdirectories, and/or change the name
if [ -n "$subdir_str" ]; then
    if [ "$subdir_str" != "${subdir_str%/}" ]; then
        # The input string is only the path, and the default name will be used
        subdir_dir="${subdir_str%/}"
    else
        str_tmp="${subdir_str%/*}"
        if [ "$str_tmp" != "" ] && [ "$str_tmp" != "$subdir_str" ]; then
            # The input string contains both the name and the path
            # Extract the path
            subdir_dir="$str_tmp"
        fi
        # Extract the name (either the full input string, or the last part of it)
        subdir_prefix="${subdir_str##*/}"
    fi
fi

[ $l_verbose -eq 1 ] && printf "\nFiles are copied into subdirectories $subdir_prefix in $subdir_dir/\n\n"

[ $l_verbose -eq 1 ] && printf "\nIncluding the following files:\n\n"

for fn in "${files_to_copy[@]}"; do
    
    [ $l_verbose -eq 1 ] && echo "$fn"
    
    for str in "${strs_to_search[@]}"; do
        
        prec="${str%${delim}*}"
        succ="${str##*${delim}}"
        
        # Consider varying amounts of white space, case insensitive
        names_to_modify+=( $( grep -ioP "\b${prec}\s*\K\w+(?=\s*${succ})" "$fn" ) )
        
    done
    
done

# Remove multiple occurrences of the same name
eval names_to_modify=( $( printf "%q\n" "${names_to_modify[@]}" | sort -u ) )

# Remove the excluded names
eval names_to_modify=( $( printf "%q\n" "${names_to_modify[@]}" "${names_to_exclude[@]}" | sort -f | uniq -iu ) )

if [ $l_verbose -eq 1 ]; then
    printf "\nRenaming the following instances:\n\n"
    for str in "${names_to_modify[@]}"; do
        echo "$str"
    done
fi

###############################################################################

# Create a new subdirectory for each chemical system, and copy the cluster dynamics files to the subdirectory

str_inc=""

for nchem in "${nchem_arr[@]}"; do
    
    subdir="${subdir_dir}/${subdir_prefix}_${nchem}"
    
    # If the directory already exists
    if [ -d "$subdir" ]; then
        
        read -p "Subdirectory $subdir already exists; do you want to write over it (y/n)? " answer
        
        if [ "$answer" != "${answer#[Yy]}" ]; then
            # Remove the contents of the previous directory
            rm -rf "$subdir/"*
        else
            printf "\nNo files copied for system $nchem.\n"
            continue
        fi
    fi
    
    mkdir -p "$subdir"
    
    # Copy the files and add the suffices and other modifications
    for fn in "${files_to_copy[@]}"; do
        
        # Insert the suffix before the file name extension
        fn_base="${fn##*/}"
        fpre="${fn_base%.*}"
        fext="${fn_base##*.}"
        fn_cp="$subdir/${fpre}_${nchem}.${fext}"
        
        cp "$fn" "$fn_cp"
        
        # Rename the routines within the file in case this is a Fortran source file
        if [ "$fext" != "${fext#[Ff]}" ]; then
            
            for str in "${names_to_modify[@]}"; do
                perl -pi -e "s/\b${str}\b/${str}_${nchem}/g" "$fn_cp"
            done
            
        # Modify the Perl call to include the suffix in the Perl output
        elif [ "$fn_base" == "$run_perl_file" ]; then
            
            for str in "${params_to_set[@]}"; do
                
                prec="${str%${delim}*}"
                succ="${str##*${delim}}"
                
                [ "$succ" == "$suffix_key" ] && succ="\"_${nchem}\""
                
                perl -pi -e "s/^[#\s]*\b${prec}\b\s*=.*$/${prec}=${succ}/g" "$fn_cp"
                
            done
            
            if [ $l_verbose -eq 1 ] && [ $nchem -eq ${nchem_arr[0]} ]; then
                printf "\nChanges in $run_perl_file compared to the original file:\n\n"
                sdiff -s "$fn" "$fn_cp"
            fi
            
        fi
        
    done
    
    # Check that the modifications have been correctly inserted (the replacements are done after all in a bit of a "black-box" manner)
    if [ $l_verbose -eq 1 ] && [ $nchem -eq ${nchem_arr[0]} ]; then
        printf "\nOccurrences of the number suffix in the copied files (messy output, but it's just to see that everything is ok):\n\n"
        grep --color "_${nchem}" "$subdir/"*
    fi
    
    # Parse also a string to be printed in a Fortran include file for automatically getting the needed modules
    str_inc+="use clusterin_${nchem}, only : cluster_dynamics_${nchem}\n"
    str_inc+="use acdc_simulation_setup_${nchem}, only : get_system_size_${nchem}\n"
    
done

printf "$str_inc" > "$subdir_dir/$inc_use_file"

printf "\nNext, you need to run $run_perl_file in each subdirectory to set up the corresponding chemical system.\n\n"
