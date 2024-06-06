
export intelDIR=/opt/intel/oneapi

export mklDIR=$intelDIR/mkl/2023.2.0/lib/intel64
echo $mklDIR

d="$(find $intelDIR -name "libmkl*" | head -1 |xargs dirname)"
export mklDIR=$d
echo $mklDIR
#find $intelDIR -name "libmkl_intel*.so"
#$intelDIR/mkl/2023.2.0/lib/intel64   libmkl_intel_lp64.so

d="$(find $intelDIR -name "libmpifort*" | head -1 |xargs dirname)"

export LIBRARY_PATH=./
export LIBRARY_PATH=$LIBRARY_PATH:$intelDIR/mpi/2021.10.0/lib/release
export LIBRARY_PATH=$LIBRARY_PATH:$intelDIR/mpi/2021.10.0/lib/
#export LIBRARY_PATH=$LIBRARY_PATH:$intelDIR/compiler/2023.2.0/linux/compiler/lib/intel64_lin/ # libimf.so

export LD_LIBRARY_PATH=./
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$mklDIR #libmkl_intel_lp64.so.2
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$intelDIR/mpi/2021.10.0/lib/  # libmpifort.so
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$intelDIR/mpi/2021.10.0/lib/release/ #libmpi.so.12
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$intelDIR/compiler/2022.1.0/linux/compiler/lib/intel64_lin #libmpi.so

export PATH=$PATH:$intelDIR/compiler/2022.1.0/linux/bin/intel64/
export PATH=$PATH:$intelDIR/mpi/2021.10.0/bin/

               FC=$intelDIR/compiler/2022.1.0/linux/bin/intel64/ifort

echo PATH=$PATH
echo LIBRARY_PATH=$LIBRARY_PATH
echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
comando="which ifort"
echo $comando; eval $comando
comando="which mpiifort"
echo $comando; eval $comando
comando="mpiifort -v"
echo $comando; eval $comando


