#!/bin/bash

set -e


# Benchmarking Batch File For Raven version evaluation
echo "benchmarking raven..."

# Location of Working directory (no end slash)
workingdir=$PWD

# reference version
ref_ver_name="ref"

# version name of NEW version
ver_name="new" #"v76"

# Location of Raven executable
ravexe=${workingdir}"/_Executables/"${ver_name}"/Raven.exe"

if [ ! -e ${ravexe} ] ; then
  echo "raven file executable "${ravexe}" doesn't exist. BENCHMARKING FAILED."
  exit 1
fi

if [ ! -e ${workingdir} ] ; then
  echo "working directory "${workingdir}" doesn't exist. BENCHMARKING FAILED."
  exit 2
fi

# create output directory for this version
if [ -e ${workingdir}"/out_"${ver_name} ] ; then
    rm -r ${workingdir}"/out_"${ver_name}
fi
mkdir ${workingdir}"/out_"${ver_name}

test_cases=(    "Alouette" "Alouette2" "York" "Irondequoit" "LOTW" "LaJoie" "Nith" "Revelstoke" "Salmon_GR4J" "Salmon_HBV" "Salmon_HMETS" "Salmon_MOHYSE" "Williston_Finlay" "York_gridded_m_daily_i_daily" "York_nongridded_m_daily_i_daily")
test_cases_rvi=("Alouette_ws" "Alouette2" "York_gridded_m_daily_i_daily" "Irondequoit" "LOWRL" "La_Joie_ws" "Nith" "Revelstoke_ws" "raven-gr4j-salmon" "raven-hbv-salmon" "raven-hmets-salmon" "raven-mohyse-salmon" "Williston_Finlay_ws" "York_gridded_m_daily_i_daily" "York_nongridded_m_daily_i_daily")

ntest_cases=$( echo "${#test_cases[@]}" )
# 

for (( icase=0 ; icase < ${ntest_cases} ; icase++ )) ; do

    test_case=${test_cases[icase]}
    rvi_file=${test_cases_rvi[icase]}
    
    echo ""
    echo ""
    echo "testing "${test_case}" ... "
    echo ""
    echo ""
    mkdir ${workingdir}"/out_"${ver_name}"/"${test_case}"/"
    cd    ${workingdir}"/_InputFiles/"${test_case}"/"

    outdir=${workingdir}"/out_"${ver_name}"/"${test_case}"/"    
    ${ravexe} ${rvi_file} -o ${outdir} > tmp.tmp 2>&1

    grep "Successful Simulation" tmp.tmp  # just to have something written to screen

    outfiles=$( \ls ${outdir}/* )

    for ff in ${outfiles} ; do

	ref_file="../../out_"${ref_ver_name}"/"${test_case}"/"$( echo $ff | rev | cut -d / -f 1 | rev )

	echo "compare to reference: "${ff}
	# echo ${ref_file}
	diff ${ff} ${ref_file}
	
    done

    if [ -e tmp.tmp ] ; then
	rm tmp.tmp
    fi

done

echo "-------------------------------------------"
echo "... BENCHMARKING DONE."
echo "-------------------------------------------"

exit 0

