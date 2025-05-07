#!/bin/bash

# Script to run V2 on the alpine-mpich Docker cluster

# Parameters
NUM_PROCS=$1
SPARSIFICATION_P=$2
PROJECT_DIR=/home/arshiq/alpine-mpich/cluster/project
DATASET=../Data/rMatGraph_J_5_100.txt
LOG_DIR=logs

# Validate inputs
if [ -z "$NUM_PROCS" ]; then
    echo "Error: Number of processes not provided."
    exit 1
fi
if [ -z "$SPARSIFICATION_P" ]; then
    echo "Error: Sparsification probability not provided."
    exit 1
fi

# Clear and create project directory
rm -rf $PROJECT_DIR/*
mkdir -p $PROJECT_DIR
mkdir -p $PROJECT_DIR/Data
mkdir -p $LOG_DIR
chmod -R 777 $LOG_DIR

# Copy files to project directory
cp v2.cpp $PROJECT_DIR/
cp Makefile $PROJECT_DIR/
cp $DATASET $PROJECT_DIR/Data/

# Change to cluster directory
cd /home/arshiq/alpine-mpich/cluster

# Start cluster
./cluster.sh $NUM_PROCS

# Wait for cluster to be ready
echo "Waiting for cluster to start..."
sleep 60

# Compile and run in master container
docker exec cluster-master-1 /bin/sh -c "cd /project && make butterfly_v2 && echo -e '2\n1\n0' | mpirun -np $NUM_PROCS ./butterfly_v2"

# Copy logs
mkdir -p /home/arshiq/Documents/GitHub/PDC_ParButterfly/src/V2/$LOG_DIR
docker cp cluster-master-1:/project/output_*.txt /home/arshiq/Documents/GitHub/PDC_ParButterfly/src/V2/$LOG_DIR/

# Stop cluster
docker compose down

echo "Cluster execution complete. Logs copied to $LOG_DIR/"