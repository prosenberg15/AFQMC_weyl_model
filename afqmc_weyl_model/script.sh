#!/bin/bash
OLDIFS=$IFS
IFS=$'\n'
export set=2           #set:1 by hand 2 TBC
export Nsite=9         #Nsite->set=1
export Nhop=36         #Nhop->set=1
export Nl1=9           #Nl(1)->set/=1 & Dimen>=1
export Nl2=9           #Nl(2)->set/=1 & Dimen>=2
export Nl3=1           #Nl(3)->set/=1 & Dimen=3
export kbound1=0.d0    #kbound(1)->set=3 & Dimen>=1
export kbound2=0.d0    #kbound(2)->set=3 & Dimen>=2
export kbound3=0.d0    #kbound(3)->set=3 & Dimen>=3
export openbcx=1     #open BCs -> 1, other BCs ->0
export openbcy=0     #open BCs -> 1, other BCs ->0
export t1=-1.d0        #t1 Hubbard hopping t1 in nearest direction.
export lamda=0.43037d0      #lamda spin orbit coupling strength
export onsitU=-5.13277d0    #onsitU Hubbard U interaction on the same site.
export Ntot=144          #Ntot the up spin and down spin particle numbers
export dtype='c'       #dtype for the determinate type: d decouple, c couple.
export Nspin1=7        #Nspin(1) means Nup
export Nspin2=7       #Nspin(2) means Ndn
export dt=0.025d0       #dt each slice of imagine time
export kcrn=1          #kcrn 1.d s;2.d c:for different kinds of release method
export bgset=1         #bgset 0. mean field 1. dynamic background walker
export pfft=0          #pfft 0.do not use fftw  1. use fftw for the code
export diagm=0         #diagm 0. get all the matrix 1. only measure diaganal matrix. pfft must be 1.
export Ntherm=15        #Ntherm number of therm
export Nmeas=60       #Nmeas number of measurement
export StepforGram=10  #StepforGram number of steps when we do modified GS
export PP=0            #PP:0 phi and aux from code, 1 aux from code read phi, 2 read aux and phi
export Nlen=1500       #Nlen for the length of the beta
export blk=50          #blk the size of the block
export meastep=25      #meastep, how often to measure.
export thermstep=500   #thermstep, therm for non-commute observables.
export Nmpi=

dir_now=`pwd`

paramfile=$1
#while read line
for line in `cat $paramfile`
do
    params=$line
    Ntot=`echo $params | awk '{print $1}'`
    Nl1=`echo $params | awk '{print $2}'`
    Nl2=`echo $params | awk '{print $3}'`
    #Nl2=$Nl1
    lamda=`echo $params | awk '{print $4}'`
    onsitU=`echo $params | awk '{print $5}'`
    TIME=`echo $params | awk '{print $6}'`

    mkdir ../result_hu_"$Nl1"_"$Nl2"_U"$onsitU"
#    echo $params | awk '{LX=$1;LY=$1;LAMBDA=$3;UHUBB=$4;TIME=$5}'

#cp  for_k_point/k_point_1 ./
#sort -k2n k_point_1 | uniq > k_point_2
#gawk '{$1=$1"d0";$2=$2"d0";print $0}' k_point_2 > k_point
#file="k_point"
#for kxky in `cat $file`
#do
 # echo $kxky| gawk '{kx=$1;ky=$2;print kx;print ky}'
#  kbound1=`echo $kxky| gawk '{print $1}'`
#  kbound2=`echo $kxky| gawk '{print $2}'`
# this is from the original script
cat >param <<!
${set}        !set:1 by hand 2 TBC
${Nsite}      !Nsite->set=1
${Nhop}       !Nhop->set=1
${Nl1}        !Nl(1)->set/=1 & Dimen>=1
${Nl2}        !Nl(2)->set/=1 & Dimen>=2
${Nl3}        !Nl(3)->set/=1 & Dimen=3
${kbound1}    !kbound(1)->set=3 & Dimen>=1
${kbound2}    !kbound(2)->set=3 & Dimen>=2
${kbound3}    !kbound(3)->set=3 & Dimen>=3
${openbcx}     !open BCs -> 1, other BCs ->0
${openbcy}     !open BCs -> 1, other BCs ->0
(${t1},0.d0)  !t1 Hubbard hopping t1 in nearest direction.
${lamda}      !lamda spin orbit coupling strength
${onsitU}     !onsitU Hubbard U interaction on the same site.
${Ntot}       !Ntot the up spin and down spin particle numbers
${dtype}      !dtype for the determinate type: d decouple, c couple.
${Nspin1}     !Nspin(1) means Nup
${Nspin2}     !Nspin(2) means Ndn
${dt}         !dt each slice of imagine time
${kcrn}       !kcrn 1.d s;2.d c for different kinds of release method
${bgset}      !bgset 0. mean field 1. dynamic background walker
${pfft}       !pfft 0.do not use fftw  1. use fftw for the code
${diagm}      !diagm 0. get all the matrix 1. only measure diaganal matrix. pfft must be 1.
${Ntherm}     !Ntherm number of therm
${Nmeas}      !Nmeas number of measurement
${StepforGram}   !StepforGram number of steps when we do modified GS
${PP}         !PP:0 phi and aux from code, 1 aux from code read phi, 2 read aux and phi
${Nlen}       !Nlen for the length of the beta
${blk}        !blk the size of the block
${meastep}    !meastep how often to measure.
${thermstep}  !thermstep therm for non-commute observables.
!

export gz=mafqmc1.8-pairing-m${Nmpi}-${set}-${Nsite}-${Nhop}-${Nl1}-${Nl2}-${Nl3}-${kbound1}-${kbound2}-${kbound3}-${openbcx}-${openbcy}-${t1}-${lamda}-${onsitU}-${Ntot}-${dtype}-${Nspin1}-${Nspin2}-${dt}-${kcrn}-${bgset}-${pfft}-${diagm}-${Ntherm}-${Nmeas}-${StepforGram}-${PP}-${Nlen}-${blk}-${meastep}-${thermstep}.tar.gz
export wkdir=mafqmc1.8-pairing-m${Nmpi}-${set}-${Nsite}-${Nhop}-${Nl1}-${Nl2}-${Nl3}-${kbound1}-${kbound2}-${kbound3}-${openbcx}-${openbcy}-${t1}-${lamda}-${onsitU}-${Ntot}-${dtype}-${Nspin1}-${Nspin2}-${dt}-${kcrn}-${bgset}-${pfft}-${diagm}-${Ntherm}-${Nmeas}-${StepforGram}-${PP}-${Nlen}-${blk}-${meastep}-${thermstep}

#tar -cf ${gz} afqmc param *.dat
tar -cf ${gz} afqmc param
#mv ${gz} ../result_"$Nl1"_"$Nl2"/
mv ${gz} ../result_hu_"$Nl1"_"$Nl2"_U"$onsitU"
#cd ../result_"$Nl1"_"$Nl2"/
cd ../result_hu_"$Nl1"_"$Nl2"_U"$onsitU"
rm -rf ${wkdir}
mkdir ${wkdir}
mv ${gz} ${wkdir}
cd ${wkdir}
tar -xf ${gz}
rm -rf ${gz}

export RUNDIR='${PBS_O_WORKDIR}'

#PBS -l nodes=1:ice:ppn=24

cat >script <<!
#!/bin/tcsh
#PBS -N afqmc
#PBS -l walltime=${TIME}:00:00
#PBS -l nodes=1:x5672:ppn=1

cd ${RUNDIR}

mvp2run ./afqmc >> out

!

qsub script
cd $dir_now
#done
IFS=$OLDIFS
echo "sbatch---done"

done < $paramfile
#nohup mpirun -np ${Nmpi} cpqmc &
#SBATCH --exclusive
#SBATCH -C new
