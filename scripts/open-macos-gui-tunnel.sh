#!/bin/bash

# http://joshuamccall.com/articles/docker.html
open -a XQuartz
export ip0=$(ifconfig en0 | grep inet | awk '$1=="inet" {print $2}')
export ip1=$(ifconfig en1 | grep inet | awk '$1=="inet" {print $2}')
export ip=$ip0$ip1
xhost + $ip

# Write to file "ip"
file=".ip"
echo "$ip" > "$file"