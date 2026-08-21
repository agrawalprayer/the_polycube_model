#!/bin/bash

# Hierarchical Pairwise Data Split Combiner and Job Manager
# [Author: Prarthana Agrawal]
# [Date: 26 Jun 2026]
# [Version: 0.3]

# Description:

    # Hierarchical pairwise combination of data splits:
    # Iteratively merges splits in pairs, doubling group size each level,
    # until all are combined into a final output.
    # Manages job submission and waits for completion at each level.
    # For tot_splits = 10
    # 1-2 3-4 5-6 7-8 9-10 (after level1 combine.py)
    # 1-4 5-8 9-10 (combining level1 output using groupcombine.py)
    # 1-8 9-10 level3
    # 1-10 level4 is the final output here -- rename to combined_*.txt

# Versions
    # 21 June 2026 v0.2: File2 shapes need be compared in serial order. We can split this into groups and run parallely.
    
    # 19 Jun 2026: Included input option to resume from a specific level. 

    # 18 May 2025: Initial version created. Basic functionality implemented.
    #    - Updates in next version: increases polling interval to 1hour.

# Usage: ./combine_ctrl.sh from the directory, or 
# Usage[UPDATED]: ./combine_ctrl.sh <level1_group_size> <workers_per_job> [start_level]
# addqueue the script using: addqueue -q long -c "controller" combine_ctrl.sh $group_size $workers_per_job $start_level

# Create combined_files directory if it doesn't exist
mkdir -p combined_files
mkdir -p outputs/combine_outputs

#---------------------------------- Import values from input.py ----------------------------------------- #
n_tiles=$(python3 -c "import input; print(input.n_tiles)")
n_sides=$(python3 -c "import input; print(input.n_sides)")
n_rules=$(python3 -c "import input; print(input.n_rules)")
tot_splits=$(python3 -c "import input; print(input.tot_splits)")
echo "Total splits: $tot_splits"
echo "System: ${n_tiles}s${n_sides}c ${n_rules}"
#-------------------------------------------------------------------------------------------------------- #

#-------------------------------------------------------------------------------------------------------- #
# Check that all required nsplit files are present
#-------------------------------------------------------------------------------------------------------- #

check_all_nsplits() {
    missing=()

    for ((i=1; i<=tot_splits; i++)); do
        if [[ ! -f "data_files/${n_tiles}s${n_sides}c_nsplit_${i}_frequency.txt" ]]; then
            missing+=("$i")
        fi
    done

    if (( ${#missing[@]} == 0 )); then
        return 0
    else
        return 1
    fi
}

echo "Waiting for all nsplit files to become available..."

while true; do
    if check_all_nsplits; then
        echo "✅ All nsplit files (1-${tot_splits}) are present. Starting combining."
        break
    fi

    echo "⏳ Still missing: ${missing[*]}"
    echo "Sleeping for 1 hour..."
    sleep 3600
done

# ========================================================================================================= #
# -------------------------------------------- Job Submission -------------------------------------------- #
# ========================================================================================================= #

# Initial parameters
base_group_size=${1:-5}
workers_per_job=${2:-20}
start_level=${3:-1}

level=$start_level
group_size=$(( base_group_size * (2 ** (start_level - 1)) ))

prev_groups=()
level1_job_ids=()
script_path="$(pwd)/run_python_job.sh"

build_groups () {
    local gsize=$1
    local groups=()
    local i end

    for ((i=1; i<=tot_splits; i+=gsize)); do
        end=$((i + gsize - 1))
        (( end > tot_splits )) && end=$tot_splits
        groups+=("$i-$end")
    done

    printf '%s\n' "${groups[@]}"
}

# If starting from Level 2 or later, reconstruct the previous level groups automatically
if (( start_level > 1 )); then
    prev_group_size=$(( base_group_size * (2 ** (start_level - 2)) ))
    mapfile -t prev_groups < <(build_groups "$prev_group_size")
fi

while true; do
    current_groups=()
    job_ids=()  # Array to store job IDs for current level
    echo -e "\n-------------------- Level $level --------------------"
    
    for ((i=1; i<=tot_splits; i+=group_size)); do
        end=$((i + group_size - 1))
        (( end > tot_splits )) && end=$tot_splits # Adjust end to not exceed tot_splits
        group_range="$i-$end" 
        current_groups+=("$group_range")

        echo "Combining #$i to #$end"
        job_name="${n_tiles}s${n_sides}c_L${level}_G${i}to${end}"
        mem=$((2 ** (level + 2)))  # 8GB at level=1, 16GB at level=2, 32GB at level=3, etc.

        # Cap memory at 20 GB
        if (( mem > 20 )); then
            mem=20
        fi

        if [[ $level -eq 1 ]]; then # Level 1 jobs: eg 1-2, 3-4, 5-6, etc.
            job_id=$(sbatch --partition=long --mem=${mem}G --job-name="$job_name" --output="outputs/combine_outputs/${job_name}-%j.out" --parsable "$script_path" combine.py input.py "$i" "$end")

            # Check if job submission was successful
            if ! [[ "$job_id" =~ ^[0-9]+$ ]]; then
                echo "⚠️ Job submission failed with $mem GB per CPU. Retrying with 30GB..."

                # Try again with higher memory (30GB)
                mem=30
                job_id=$(sbatch --partition=long --mem=${mem}G --job-name="$job_name" --output="outputs/combine_outputs/${job_name}-%j.out" --parsable "$script_path" combine.py input.py "$i" "$end")

                if ! [[ "$job_id" =~ ^[0-9]+$ ]]; then
                    echo "❌ Job submission failed again at 30GB. Exiting."
                    echo "❌ Failed at Level $level for group $i-$end (job name: $job_name)."
                    exit 1  #Exit
                fi
            fi

            level1_job_ids+=("$job_id")

        else
            tuple_parts=()
            for pg in "${prev_groups[@]}"; do
                IFS='-' read -r p_start p_end <<< "$pg"
                if (( p_start >= i && p_end <= end )); then
                    tuple_parts+=("${p_start}_${p_end}")
                fi
            done

            tuple_string=$(IFS=.; echo "${tuple_parts[*]}")  # Join with period as delimiter
            tuple_string="${tuple_string// /}"  # Remove any whitespace

            job_id=$(sbatch --partition=long --cpus-per-task="$workers_per_job" --mem=${mem}G --job-name="$job_name" --output="outputs/combine_outputs/${job_name}-%j.out" --parsable "$script_path" groupcombine.py input.py "$tuple_string")
            
            echo "Job ID: $job_id"

            # Check if job submission was successful
            if ! [[ "$job_id" =~ ^[0-9]+$ ]]; then
                echo "⚠️ Job submission failed with $mem GB per CPU. Retrying with 30GB..."

                # Try again with higher memory (30GB)
                mem=30
                job_id=$(sbatch --partition=long --cpus-per-task="$workers_per_job" --mem=${mem}G --job-name="$job_name" --output="outputs/combine_outputs/${job_name}-%j.out" --parsable "$script_path" groupcombine.py input.py "$tuple_string")
                echo "Job ID: $job_id"

                if ! [[ "$job_id" =~ ^[0-9]+$ ]]; then
                    echo "❌ Job submission failed again at 30GB. Exiting."
                    echo "❌ Failed at Level $level for $tuple_string (job id: $job_id)."
                    exit 1  #Exit
                fi
            fi

            job_ids+=("$job_id")
        fi
    done

    # ======================================================================================================= #
    # -------------------------------------------- Job Monitoring ------------------------------------------- #
    # ======================================================================================================= #

    # ----------------------------------------------------------------------------------- #
    # Wait for Level 1 jobs to finish before continuing
    #----------------------------------------------------------------------------------- #
    if [[ $level -eq 1 ]]; then
        echo "Waiting for all Level 1 jobs to finish..."

        declare -A warned_missing_log=()  # Track jobs with missing logs to avoid spamming
        check_count=0
        while true; do
            sleep 1200  # Polling interval
            ((check_count++)) # this is to print only after check_count times, to avoid spamming the console too much

            active_jobs=()
            for jid in "${level1_job_ids[@]}"; do
                # Check if job is still in queue
                if squeue -j "$jid" > /dev/null 2>&1 && squeue -j "$jid" | grep -q "$jid"; then
                    active_jobs+=("$jid")
                else
                    read -r state exitcode runtime <<< $(sacct -j "$jid" --format=State,ExitCode,Elapsed --noheader 2>/dev/null | head -n1 | awk '{gsub(/[ \t]+/, " "); print $1, $2, $3}')
                    exitcode=${exitcode%%:*}

                    if [[ -z "$state" || "$state" == "Unknown" ]]; then
                        job_log_file="python-${jid}.out"
                        if [[ -f "$job_log_file" ]]; then
                            if grep -qiE "Traceback|error|exception|oom_kill|OOM" "$job_log_file"; then
                                echo "❌ Job $jid likely failed based on log contents."
                                echo "🔍 See $job_log_file for details. (Level $level, Job ID: $jid)"
                                exit 1
                            fi
                        else
                            # Only print warning once per missing log file to avoid spam
                            if [[ -z "${warned_missing_log[$jid]}" ]]; then
                                echo "⚠️ Log file not found for Job $jid. Proceeding cautiously."
                                warned_missing_log[$jid]=1
                            fi
                        fi
                    elif [[ "$state" == "FAILED" || "$state" == "CANCELLED" || "$exitcode" -ne 0 ]]; then
                        echo "❌ Job $jid failed with state=$state, exit code=$exitcode, runtime=${runtime:-N/A} (Level $level)"
                        exit 1
                    else
                        echo "✅ Job $jid completed successfully. Runtime: ${runtime:-N/A}"
                    fi
                fi
            done

            if [[ ${#active_jobs[@]} -eq 0 ]]; then
                echo "✅ All Level 1 jobs are completed. Proceeding to Level 2."
                break
            #else
            #    echo "⏳ Still waiting on Level 1 jobs: ${active_jobs[*]}"
            fi

            if (( check_count % 5 == 0 )); then
                echo "⏳ Still waiting on Level 1: ${#active_jobs[@]} jobs remaining"
            fi
        done
    fi

    # ---------------------------------------------------------------------------------- #
    # Wait for jobs in Level 2+ (if any)
    # ---------------------------------------------------------------------------------- #
    if [[ $level -gt 1 ]]; then
        echo "Waiting for jobs in Level $level to finish: ${job_ids[*]}"

        declare -A warned_missing_log_level=()  # Same for higher levels
        check_count=0
        while true; do
            sleep 1200  # Polling interval
            ((check_count++)) # this is to print only after check_count times, to avoid spamming the console too much
            active_jobs=()

            for jid in "${job_ids[@]}"; do
                if squeue -j "$jid" > /dev/null 2>&1 && squeue -j "$jid" | grep -q "$jid"; then
                    active_jobs+=("$jid")
                else
                    read -r state exitcode runtime <<< $(sacct -j "$jid" --format=State,ExitCode,Elapsed --noheader 2>/dev/null | head -n1 | awk '{gsub(/[ \t]+/, " "); print $1, $2, $3}')
                    exitcode=${exitcode%%:*}

                    if [[ -z "$state" || "$state" == "Unknown" ]]; then
                        job_log_file="python-${jid}.out"
                        if [[ -f "$job_log_file" ]]; then
                            if grep -qiE "Traceback|error|exception|oom_kill|OOM" "$job_log_file"; then
                                echo "❌ Job $jid likely failed based on log contents."
                                echo "🔍 See $job_log_file for details. (Level $level, Job ID: $jid)"
                                exit 1
                            fi
                        else
                            # Only warn once per missing log to reduce spam
                            if [[ -z "${warned_missing_log_level[$jid]}" ]]; then
                                echo "⚠️ Log file not found for Job $jid. Proceeding cautiously."
                                warned_missing_log_level[$jid]=1
                            fi
                        fi
                    elif [[ "$state" == "FAILED" || "$state" == "CANCELLED" || "$exitcode" -ne 0 ]]; then
                        echo "❌ Job $jid failed with state=$state, exit code=$exitcode, runtime=${runtime:-N/A} (Level $level)"
                        exit 1
                    else
                        echo "✅ Job $jid completed successfully. Runtime: ${runtime:-N/A}"
                    fi
                fi
            done

            if [[ ${#active_jobs[@]} -eq 0 ]]; then
                echo "✅ All jobs for Level $level completed."
                break
            #else
            #    echo "⏳ Still waiting on jobs in Level $level: ${active_jobs[*]}"
            fi

            if (( check_count % 5 == 0 )); then
                echo "⏳ Still waiting on Level $level: ${#active_jobs[@]} jobs remaining"
            fi
        done
    fi

    # Stop condition
    if (( group_size > tot_splits )); then
        break
    fi
    prev_groups=("${current_groups[@]}")
    group_size=$((group_size*2))
    ((level++))
done

# Rename final grouped files
echo "Renaming final grouped files to combined_*"
combined_files="combined_files"

for file in "$combined_files/${n_tiles}s${n_sides}c_combined_1to${tot_splits}nsplits"*.txt; do
    [[ -f "$file" ]] && mv "$file" "${file/_1to${tot_splits}nsplits/}"
done
