python ../ConvertCellModel.py "$1.cellml" -tOdeint
mv "$1.cpp" "$1.cu"
nvcc -o "$1.out"  "$1.cu" -I /home/ug13ag2/odeint/include
