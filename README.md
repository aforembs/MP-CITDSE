# 2e-DMX
2 electron dipole matrices from 1 electron partial waves

## Requirements
A C++17 or greater compatible compiler
tbb  - library for C++17 multithreading (apt-get, pacman, etc...)
hdf5 - data storage library (apt-get, pacman, etc...)

## Installation 
To compile do 
`cd src`
`mkdir obj`
`make`

To run do
`./dmx_test`

This produces he0idx.h5 , he1idx.h5, ... index files 
and he2\_01v.h5, he2\_12v.he ... 2e dmx files
