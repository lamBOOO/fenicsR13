#!/bin/bash

# http://joshuamccall.com/articles/docker.html
open -a XQuartz
export ip=$(ifconfig en0 | grep inet | awk '$1=="inet" {print $2}')
export TEST=hi
xhost + $ip