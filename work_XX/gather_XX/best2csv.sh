#!/bin/sh

grep best *.gather >& energy.XX_mm.bcc_allgather.best

sed 's/  / /g' < energy.XX_mm.bcc_allgather.best > ${1}
sed -i 's/  / /g' ${1}
sed -i 's/  / /g' ${1}
sed -i 's/  / /g' ${1}
sed -i 's/ /,/g' ${1}
