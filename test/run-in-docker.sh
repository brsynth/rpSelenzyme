#!/bin/bash


cd ../docker
docker-compose \
  --compatibility run \
  --rm \
  -v $PWD/../test:/home/test \
  -w /home/test --entrypoint="" \
  tool \
    sh -c "./run.sh"
cd -
