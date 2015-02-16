#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
set -e
trap 'kill $(jobs -p)' SIGINT SIGTERM EXIT
export TRAVIS=1
export DOCKER_HOST="unix:///var/run/docker.sock"

docker -d &
sleep 2
chmod +rw /var/run/docker.sock

docker version

docker pull satta/annot-nf
ls -Al
./nextflow -c loc_docker.config -c params_default.config run annot.nf -with-docker satta/annot-nf

killall -9 docker