#!/bin/bash
# Example: "repack_hdf5_files.sh ." will repack all files matching tracer_file*.hdf5 in the cwd

# for filename in $1/tracer_file*.hdf5; do
# 	echo "Repacking file: $filename"
# 	h5repack "$filename" "tmp.hdf5"
# 	rm "$filename"
# 	mv "tmp.hdf5" "$filename"
# done

for i in {seq -f "%03g" $2 $3}; do
	filename = $1/tracer_file_$i.hdf5
	echo "Repacking file: $filename"
	h5repack "-f GZIP=9 $filename tmp.hdf5"
	rm "$filename"
	mv "tmp.hdf5" "$filename"
done
