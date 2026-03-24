#!/bin/bash
#SBATCH --output=python-%j.out # Standard output file (%j = job ID)

# Record start time
start_time=$(date +%s)

# Optional: print job info
echo "Running job on $(hostname) at $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "------------------------------------------------------------"

# Activate your virtual environment and run the Python script passed as an argument
/mnt/users/pagrawal/polycube/.venv/bin/python "$@"

# Check the exit code of the Python script
exit_code=$?

echo "Python script exited with code $exit_code"

# Record end time and calculate elapsed time
end_time=$(date +%s)
elapsed=$(( end_time - start_time ))

# Format elapsed time in hh:mm:ss
printf "Total run time: %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
echo "------------------------------------------------------------"

exit $exit_code
