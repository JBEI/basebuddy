#!/bin/bash

echo "Running BaseBuddy with Docker locally."

if [[ ! -n $1 ]];
then 
    PORT=8080
else
    PORT=$1
fi

echo "Connect using port $PORT.\n\n"

docker build . -t basebuddy
docker run -v "${pwd}:/app:ro" -p "$PORT:8080" -i basebuddy