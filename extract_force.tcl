# Open the input log file for reading
set input_file [open "9A_002M_PME_D12_SMD_4ns.6295880-out" "r"]

# Open the output file to write data
set output_file [open "ft_2.dat" "w"]

# Gather input from the user
puts "Enter a value for n_x:"
set nx [gets stdin]
puts "Enter a value for n_y:"
set ny [gets stdin]
puts "Enter a value for n_z:"
set nz [gets stdin]

# Loop through the log file
while {[gets $input_file line] != -1} {
    # Check if a line contains SMD output
    if {[lindex $line 0] eq "SMD"} {
        set step [lindex $line 1]
        # Use expr for arithmetic operations directly
        set dot_product [expr {$nx * [lindex $line 5] + $ny * [lindex $line 6] + $nz * [lindex $line 7]}]
        puts $output_file "$step $dot_product"
    }
}

# Close the input and output files
close $input_file
close $output_file
