#!/bin/bash 

datafile=input.dat
sub=job.sge
for U1 in 0.0 #1.0 2.0 3.0 4.0 5.0 6.0 #8.0 10.0  12.0 14.0 16.0 18.0 20.0  #2.0 #  8.0 10.0 12.0 14.0 16.0 18.0 20.0 
do
mkdir -p U_$U1
cd U_$U1

for tp in 0.0 #0.2 0.4 0.6 0.8 1.0 
do 
mkdir -p tp_$tp
cd tp_$tp

for seed in 87959 #98352 79869 69857 96387
do
 
mkdir -p seed_$seed

cd seed_$seed

cp ../../../ionic_hub_tca_mf_mc.f90 .
cp ../../../mtfort90.f .


#!******************************************************************************

echo   "8                          !d "             > $datafile
echo   "4                          !ds "            >> $datafile
echo   "-1.0d0                     !t1     "        >> $datafile
echo   $tp"0d0                      !t2     "        >> $datafile
echo   "14                         !temp_max "      >> $datafile
echo   "600                        !MCSW "          >> $datafile
echo   "5                          !intrvl   "      >> $datafile
echo   $U1"d0                      !U1 "             >> $datafile  
echo   "1.0d0                      !filling "       >> $datafile
echo   "0.1                       !gama "          >> $datafile 
echo   "0.1                       !gama_m "        >> $datafile 
echo   "$seed                      !seed "          >> $datafile
echo   "1.0                        !strnth  "       >> $datafile
#!******************************************************************************

gfortran -o run.x  mtfort90.f ionic_hub_tca_mf_mc.f90 -llapack -lblas
./run.x <input.dat> main_out &
#********************************************************************************
echo "#PBS  -N   U_"$U1"_tp_"$tp                   > $sub
echo "PBS   -l nodes=1:ppn=1"                      >> $sub         
echo "PBS   -koe"                                  >> $sub                
echo "cd"    \$PBS_O_WORKDIR                       >> $sub
echo "date"                                        >> $sub
echo "./run.x <input.dat |tee out "                >> $sub
echo "date"                                        >> $sub
#******************************************************************************
chmod 777 job.sge
#qsub  job.sge
cd ..
done

cd ..
done 


cd ..
done
