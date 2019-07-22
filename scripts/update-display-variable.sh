#!/bin/bash

file=".ip"
ip=$(cat "$file")

export DISPLAY="$ip":0