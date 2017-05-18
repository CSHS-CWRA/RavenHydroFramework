#!/bin/bash

set -e


# Benchmarking Batch File For Raven version evaluation
echo "benchmarking raven..."

# Location of Working directory (no end slash)
workingdir=$PWD

# version name 
ver_name="v60"

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

test_cases="York York2 Alouette Alouette2 La_Joie Revelstoke Williston_Finlay Irondequoit Nith"

for test_case in ${test_cases} ; do
    echo ""
    echo ""
    echo "testing "${test_case}" ... "
    echo ""
    echo ""
    mkdir ${workingdir}"/out_"${ver_name}"/"${test_case}"/"
    cd    ${workingdir}"/_InputFiles/"${test_case}"/"
    if [[ ${test_case} == 'Alouette' || ${test_case} == 'La_Joie' || ${test_case} == 'Revelstoke' || ${test_case} == 'Williston_Finlay' ]] ; then
	${ravexe} ${test_case}"_ws" -o ${workingdir}"/out_"${ver_name}"/"${test_case}"/"
	
    else
	if [[ ${test_case} == 'York' ]] ; then
	    ${ravexe} ${test_case}"_gridded_m_daily_i_daily" -o ${workingdir}"/out_"${ver_name}"/"${test_case}"/"	
	else
	    if [[ ${test_case} == 'York2' ]] ; then
		${ravexe} ${test_case}"_gridded_m_subdaily_i_subdaily" -o ${workingdir}"/out_"${ver_name}"/"${test_case}"/"	 
	    else
		${ravexe} ${test_case} -o ${workingdir}"/out_"${ver_name}"/"${test_case}"/"
	    fi
	fi
    fi
done

echo "-------------------------------------------"
echo "... BENCHMARKING DONE."
echo "-------------------------------------------"

exit 0

