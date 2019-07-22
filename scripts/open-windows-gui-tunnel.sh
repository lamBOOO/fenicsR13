#!/bin/bash

/c/Program\ Files/VcXsrv/xlaunch.exe -run etc/config.xlaunch
export ip=$(hostname)

# Write to file "ip"
file=".ip"
echo "$ip" > "$file"