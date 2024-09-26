## for running 'run_decode_barcode.py' on multiple fovs

# run 'run_decode_barcode.py' on one fov
#python run_decode_barcode.py 20 'fov-20'
#echo "Processing completed."

# run 'run_decode_barcode.py' over all fovs
# arguments:    first argument specifies the fov by the order in which it is read in 'run_decode_barcode.py'
#                   desired fovs are specified in the range of the 'for' loop
#               second argument names the subfolder the results are saved to
for i in {70..189}
do 
    python run_decode_barcode.py $i 'all-fovs'
done
echo "All fovs processed."

# takes about 60-70 min for all 190 fovs (20 sec each)

## stopped at 71 (70)