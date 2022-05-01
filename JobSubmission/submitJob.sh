if [ $# -lt 1 ]; then
  echo -e "\nJob script for submission of analysis jobs on KRONOS.\n"
  echo -e "USAGE: sbatch --array=1-n --$PWD/Submit.sh <input_file_list>  <prefix>\n"
  exit 1
fi

input_list=$1
prefix=$2

readarray list < $input_list
file=${list[$SLURM_ARRAY_TASK_ID-1]}
input=$file
name=$prefix"_"$SLURM_ARRAY_TASK_ID
./analysisScript.sh $input $name

echo "Job submitted to batch farm!"
