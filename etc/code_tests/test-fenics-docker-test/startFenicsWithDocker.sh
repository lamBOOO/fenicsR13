#!/bin/bash

# share webview and vurrent directory
docker run -ti -p 127.0.0.1:8000:8000 -v $(pwd):/home/fenics/shared quay.io/fenicsproject/stable:latest