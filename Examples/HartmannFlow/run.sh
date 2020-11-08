#!/bin/bash
FILE=Flow
FILEOUT=${FILE}
CHECKFILESDIR=checkfiles

rm *.chk
rm *.vtu
rm *.fld
rm -R checkfiles

/usr/bin/time -v WeakMHDSolver.out ${FILE}.xml

mkdir ${CHECKFILESDIR}
mv *.chk ${CHECKFILESDIR}

cd ${CHECKFILESDIR}
for filechk in ${FILE}_*.chk 
do
filechk=$(echo $filechk | sed 's:\(.*\)\.chk:\1:')
filechkout=$(echo $filechk | sed s:${FILE}:${FILEOUT}:)
FieldConvert ../${FILE}.xml ${filechk}.chk ${filechkout}.vtu
done
cd ..
FieldConvert ${FILE}.xml ${FILE}.fld ${FILEOUT}.vtu


# paraview --data=${FILEOUT}.vtu
