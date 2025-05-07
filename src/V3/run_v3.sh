#!/bin/bash

# Script to run V3 on the alpine-mpich Docker cluster

# Parameters
PROJECT_DIR=/home/arshiq/alpine-mpich/cluster/project
DATASET=../Data/rMatGraph_J_5_100.txt
LOG_DIR=logs

# Clear and create project directory
rm -rf $PROJECT_DIR/*
mkdir -p $PROJECT_DIR
mkdir -p $PROJECT_DIR/Data
mkdir -p $LOG_DIR
chmod -R 777 $LOG_DIR

# Copy files to project directory
cp v3.cpp $PROJECT_DIR/
cp Makefile $PROJECT_DIR/
cp $DATASET $PROJECT_DIR/Data/

# Change to cluster directory
cd /home/arshiq/alpine-mpich/cluster

# Start cluster
./cluster.sh 1

# Wait for cluster to be ready
echo "Waiting for cluster to start..."
sleep 60

# Compile and run in master container
docker exec cluster-master-1 /bin/sh -c "cd /project && make butterfly_v3 && echo -e '2\n1\n0' | ./butterfly_v3"

# Copy logs
mkdir -p /home/arshiq/Documents/GitHub/PDC_ParButterfly/src/V3/$LOG_DIR
docker cp cluster-master-1:/project/output_*.txt /home/arshiq/Documents/GitHub/PDC_ParButterfly/src/V3/$LOG_DIR/

# Stop cluster
docker compose down

echo "Cluster execution complete. Logs copied to $LOG_DIR/"