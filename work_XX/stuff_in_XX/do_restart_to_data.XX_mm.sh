#!/bin/tcsh


sed "s/Q/${1}/g" < XX_mm.restart_to_data.in > ${1}.restart_to_data.in

../source/lammps_20061129_modified/lmp_serial -in ${1}.restart_to_data.in -log  ${1}.restart_to_data.log >& ${1}.restart_to_data.out
