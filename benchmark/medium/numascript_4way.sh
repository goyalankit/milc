#!/bin/bash

# First determine "wayness" of PE
PE=4way
myway=`echo $PE | sed s/way//`

# Determine local compute node rank number
# Note, bug in MVAPICH2 uses PMI_RANK rather than MPI_RANK

if [ "x$PMI_RANK" == "x" ]; then
    local_rank=$(( $MPIRUN_RANK % $myway ))
    PMI_RANK=$MPIRUN_RANK
else
    local_rank=$(( $PMI_RANK % $myway ))
fi

# Based on "wayness" determine socket layout on local node
# if less than 4-way, offset to skip socket 0
if [ $myway -lt 4 ]; then
    numnode=$(( $local_rank + 1 ))
# if 4-way to 12-way, spread processes equally on sockets
elif [ $myway -lt 13 ]; then
    numnode=$(( $local_rank / ( $myway / 4 ) ))
# if 16-way, spread processes equally on sockets
elif [ $myway -eq 16 ]; then
    numnode=$(( $local_rank / ( $myway / 4 ) ))
# Offset to not use 4 processes on socket 0
else
    numnode=$(( ($local_rank + 1) / 4 ))
fi
echo "Running $PMI_RANK on socket $numnode"

export MV2_USE_AFFINITY=0
export MV2_ENABLE_AFFINITY=0
export VIADEV_USE_AFFINITY=0
export VIADEV_ENABLE_AFFINITY=0

# numactl -c $numnode -m $numnode numactl --show

numactl -c $numnode -m $numnode $*
