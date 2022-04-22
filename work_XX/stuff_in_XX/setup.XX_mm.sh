#!/bin/sh

if [ -s ./incoming/gb_pbv.${1} ]; then

  ln -s ../source/build_gb_diwo/gb_energy_cc_diwo .
  ln -s ../source/build_gb_diwo/build_gb_diwo .
  ln -s ../source/build_gb_diwo/gb_energy_gather .

  mkdir run_XX_mm.${1}
  cd run_XX_mm.${1}
  ln -s ../build_gb_diwo .
  sed "s/Q/XX_mm.shell.${1}/g" < ../energy.XX_mm.Q.scr > energy.XX_mm.${1}.scr
  ln -s ../energy_XX_mm.min.in .
  sed "s/Q/${1}/" < ../energy.XX_mm.Q.in > energy.XX_mm.${1}.in
  mv ../incoming/gb_pbv.${1} .

  ../gb_energy_cc_diwo energy.XX_mm.${1}.in < gb_pbv.${1} >& energy.XX_mm.${1}.out

else
  echo "GB file not found"
fi

