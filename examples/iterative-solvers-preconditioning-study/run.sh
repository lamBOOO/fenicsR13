#!/bin/bash

for p in 0 0 0 0 # TODO
do
  fenicsR13 "input$p".yml > "output$p".log
done
