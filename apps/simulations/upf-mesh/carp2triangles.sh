#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 <mesh_prefix>"
    exit
fi

mesh_prefix=$1

#
# Create .node file
#
original_file=${mesh_prefix}.pts

if [ ! -f "$original_file" ]; then
    echo "File $original_file does not exist!"
    exit
fi

new_file=${mesh_prefix}.node

# Header of the file
num_nodes=`head -n 1 $original_file`
echo "$num_nodes 3 0 0" > $new_file

# Coordinates in um should be converted to cm. Add index starting from
awk '{if (NF == 3) printf("%d %f %f %f\n", NR-2, $1/10000, $2/10000, $3/10000)}' $original_file >> $new_file


#
# Create .ele file
#
original_file=${mesh_prefix}.tetras

if [ ! -f "$original_file" ]; then
    echo "File $original_file does not exist!"
    exit
fi

new_file=${mesh_prefix}.ele

# Header of the file
num_nodes=`head -n 1 $original_file`
echo "$num_nodes 4 0" > $new_file

# Add index starting from 
awk '{if (NF == 4) printf("%d %d %d %d %d\n", NR-2, $1, $2, $3, $4)}' $original_file >> $new_file


#
# Create empty .face file
#
new_file=${mesh_prefix}.face

echo "0" > $new_file

#
# Create .fibres file
#
original_file=${mesh_prefix}.ttlon

if [ ! -f "$original_file" ]; then
    echo "File $original_file does not exist!"
    exit
fi

new_file=${mesh_prefix}.fibres

# Header of the file
num_elems=`cat $original_file | wc -l`
echo "$num_elems" > $new_file

# Fibre definition
cat $original_file >> $new_file
