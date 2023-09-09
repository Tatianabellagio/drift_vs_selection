# Paths to the input text files
tree_file="simulation_selection_drift.trees"
optima_file="optimas.txt"
variances_file="variances.txt"
# Path to the SLiM script
slim_script="arabidopsis_evolve_treeseq_d_vs_selection.slim"

echo "tree file: $tree_file" 

# Iterate through variances
while IFS= read -r variance; do
    # Iterate through optima
    while IFS= read -r optima; do
        # Construct the output file name
        
        # Construct the SLiM command
        slim_command="slim -d \"tree='$tree_file'\" -d \"optima_value='$optima'\" -d \"variance='$variance'\" \"$slim_script\""

        # Echo the command
        echo "Running: $slim_command"

        # Uncomment the following line to execute the SLiM command
        eval "$slim_command"
    done < "$optima_file"
done < "$variances_file"