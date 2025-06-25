#!/bin/bash

localport=$1

jupyter lab --no-browser &

while ! jupyter lab list | grep -q "localhost"; do
    sleep 1
done

url=$(jupyter lab list | grep -oP "http://\K[^ ]+" | tail -n 1)
port=$(echo ${url} | grep -oP "localhost:\K\d+")
statichostname=$(hostnamectl | grep -oP "Static hostname: \K.*")

if [ -z "$localport" ]; then
    localport=$port
fi

url=$(echo ${url} | sed "s/${port}/${localport}/")

echo "---------------------------------------------------
SSH in another terminal with the following command:

ssh -L ${localport}:localhost:${port} ${USER}@${statichostname}

Then, open the following link in your browser:

${url}

To stop, run 'jupyter lab stop ${port}' in the
session where you started the server.
---------------------------------------------------"