# Paths to the input text files
tree_file="../simulation_selection_drift.trees"
# Path to the SLiM script
slim_script="tree_after_cloning.slim"

# Construct the SLiM command
slim_command="slim -d \"tree='$tree_file'\" \"$slim_script\""
    
# Echo the command
echo "Running: $slim_command"
    
# Uncomment the following line to execute the SLiM command
eval "$slim_command"