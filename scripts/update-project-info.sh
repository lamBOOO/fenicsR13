#!/bin/bash

echo "Checking project's information"

commits=$(git rev-list --all --count)
linesofcode=$(find ./src -name '*.py' | xargs wc -l | tail -n1 | egrep -o "[0-9]*")
latest_tag=$(git describe --tags `git rev-list --tags --max-count=1`)

echo "{\"commits\":\"${commits}\",\"linesofcode\":\"${linesofcode}\",\"latest_tag\":\"${latest_tag}\"}" > /tmp/info.json

diff /tmp/info.json info.json

if [ $? -ne 0 ]
then
  echo "Status of project DID change (./status.json)."
  cat /tmp/info.json > info.json
  exit 1
else
  echo "Status of project did not change."
fi
