# Barnes Hut

A Barnes-hut simulator written in Regent.

## Requirements

1. A working build of Regent with HDF5 support enabled
1. (Optional) Gasnet for multi-node simulations

## Getting started

1. Convert input bodies csv to hdf5 format via hdf5_converter.rg (see sample datasets in input folder for desired format):
    - <path to regent installation>/regent.py hdf5_converter.rg -i <path_to_input_csv> -o <desired_path_to_output_hdf5_file> -n <number_of_bodies_in_input_file>
1. Run Barnes-Hut simulation:
    - Using the regent interpreter directly: <path to regent installation>/regent.py barnes_hut.rg -i <path_to_input_hdf_file> -n <number_of_bodies_in_input_file> <other options>
    - Compiling first: SAVEOBJ=1 <path to regent installation>/regent.py barnes_hut.rg, and then call the built binary with the same options as above
