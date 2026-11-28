#! /bin/bash -x

# Arguments
MPIEXEC=${1}
BINDIR=${2}
PPNAME=${3}
INITNAME=${4}
BINNAME=${5}
POSTNAME=${6}
PPCONF=${7}
INITCONF=${8}
RUNCONF=${9}
POSTCONF=${10}
PROCS=${11}
PROCS_POST=${12}
eval DATPARAM=(`echo ${13} | tr -s '[' '"' | tr -s ']' '"'`)
eval DATDISTS=(`echo ${14} | tr -s '[' '"' | tr -s ']' '"'`)

PROCLIST=(`echo ${PROCS} | tr -s ',' ' '`)
TPROC=${PROCLIST[0]}
for n in ${PROCLIST[@]}
do
   (( n > TPROC )) && TPROC=${n}
done
PROCLIST_POST=(`echo ${PROCS_POST} | tr -s ',' ' '`)

FILES_LLIO=""
if [ ! ${PPCONF} = "NONE" ]; then
   CONFLIST=(`echo ${PPCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   FILES_LLIO=`echo -e ${FILES_LLIO} ${BINDIR}/${PPNAME}`
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      RUN_PP=`echo -e "${RUN_PP}\n"${MPIEXEC} ${PROCLIST[i]} ${BINDIR}/${PPNAME} ${CONFLIST[i]} "|| exit 1"`
      FILES_LLIO=`echo -e ${FILES_LLIO} ${CONFLIST[i]}`
   done
fi

if [ ! ${INITCONF} = "NONE" ]; then
   CONFLIST=(`echo ${INITCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   FILES_LLIO=`echo -e ${FILES_LLIO} ${BINDIR}/${INITNAME}`
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      RUN_INIT=`echo -e "${RUN_INIT}\n"${MPIEXEC} ${PROCLIST[i]} ${BINDIR}/${INITNAME} ${CONFLIST[i]} "|| exit 1"`
      FILES_LLIO=`echo -e ${FILES_LLIO} ${CONFLIST[i]}`
   done
fi

if [ ! ${RUNCONF} = "NONE" ]; then
   CONFLIST=(`echo ${RUNCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   FILES_LLIO=`echo -e ${FILES_LLIO} ${BINDIR}/${BINNAME}`
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      RUN_MAIN=`echo -e "${RUN_MAIN}\n"${MPIEXEC} ${PROCLIST[i]} ${BINDIR}/${BINNAME} ${CONFLIST[i]} "|| exit 1"`
      FILES_LLIO=`echo -e ${FILES_LLIO} ${CONFLIST[i]}`
   done
elif [ ! ${BINNAME} = "NONE" ]; then
   FILES_LLIO=`echo -e ${FILES_LLIO} ${BINDIR}/${BINNAME}`
   RUN_MAIN=`echo -e "${RUN_MAIN}\n"${MPIEXEC} ${TPROC} ${BINDIR}/${BINNAME}"|| exit 1"`
fi

if [ ! ${POSTCONF} = "NONE" ]; then
   CONFLIST=(`echo ${POSTCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   FILES_LLIO=`echo -e ${FILES_LLIO} ${BINDIR}/${POSTNAME}`
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      RUN_POST=`echo -e "${RUN_POST}\n"${MPIEXEC} ${PROCLIST_POST[i]} ${BINDIR}/${POSTNAME} ${CONFLIST[i]} "|| exit 1"`
      FILES_LLIO=`echo -e ${FILES_LLIO} ${CONFLIST[i]}`
   done
fi

array=( `echo ${TPROC} | tr -s 'x' ' '`)
x=${array[0]}
y=${array[1]:-1}
let xy="${x} * ${y}"


cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# ------ For Miyabi-g
#
################################################################################
#PBS -N scale-dg
#PBS -W group_list=${Miyabi_g_GROUP}
#PBS -l select=1:mpiprocs=1:ompthreads=1
#PBS -l walltime=00:10:00
#PBS -j oe

export OMP_NUM_THREADS=1

${ADDITIONAL_CONF}

module purge
module load nvidia
module load nv-hpcx
module load cuda
module load netcdf
module load netcdf-fortran
module load hdf5

ulimit -s unlimited

EOF1

# link to file or directory
ndata=${#DATPARAM[@]}

if [ ${ndata} -gt 0 ]; then
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"

      pair=(${DATPARAM[$i]})

      src=${pair[0]}
      dst=${pair[1]}
      if [ "${dst}" = "" ]; then
         dst=${pair[0]}
      fi

      if [ -f ${src} ]; then
         echo "ln -svf ${src} ./${dst}" >> ./run.sh
      elif [ -d ${src} ]; then
         echo "rm -f          ./${dst}" >> ./run.sh
         echo "ln -svf ${src} ./${dst}" >> ./run.sh
      else
         echo "datafile does not found! : ${src}"
         exit 1
      fi
   done
fi

# link to distributed file
ndata=${#DATDISTS[@]}

if [ ${ndata} -gt 0 ]; then
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"

      triple=(${DATDISTS[$i]})

      for np in `seq 1 ${triple[0]}`
      do
         let "ip = ${np} - 1"
         PE=`printf %06d ${ip}`

         src=${triple[1]}.pe${PE}.nc
         dst=${triple[2]}.pe${PE}.nc

         if [ -f ${src} ]; then
            echo "ln -svf ${src} ./${dst}" >> ./run.sh
         else
            echo "datafile does not found! : ${src}"
            exit 1
         fi
      done
   done
fi

cat << EOF2 >> ./run.sh

# run

${RUN_PP}
${RUN_INIT}
${RUN_MAIN}
${RUN_POST}

################################################################################
EOF2
