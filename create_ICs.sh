long_x=1.0
long_y=1.0
long_z=1.0
box_size=1000000
num_tracer_x=130
num_tracer_y=1
num_tracer_z=1
output_dir="../../Scripts"
python3 create_tracer_ics.py --box_size $box_size --long_x $long_x --long_y $long_y --long_z $long_z --num_tracer_x $num_tracer_x --num_tracer_y $num_tracer_y --num_tracer_z $num_tracer_z --output_dir $output_dir
