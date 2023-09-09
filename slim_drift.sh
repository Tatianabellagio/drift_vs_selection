# Paths to the input text files
tree_file="simulation_selection_drift.trees"
# Path to the SLiM script
slim_script="arabidopsis_evolve_treeseq_drift_vs_s.slim"

#output_directory="drift_runs"

# Create the output directory if it doesn't exist

# Loop through drift_run values from 1 to 10
for drift_run in {1..10}; do
    
    # Construct the SLiM command
    slim_command="slim -d \"tree='$tree_file'\" -d \"drift_run='$drift_run'\" \"$slim_script\" "
    
    # Echo the command
    echo "Running: $slim_command"
    
    # Uncomment the following line to execute the SLiM command
    eval "$slim_command"
done