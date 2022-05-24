



def write_output(filename,dirname,submit=False):
    output = ""
    output += "#!/bin/bash"
    output += "\n"
    output += "# Set Job Name"
    output += "#SBATCH -J job"
    output += "\n"
    output += "# Set file to capture standard out and standard error and append the jobID (%j)"
    output += "#SBATCH -o job.out.%j"
    output += "\n"
    output += "#SBATCH --partition=normal"
    
    output += "#Set the number of nodes"
    output += "#SBATCH -N 1"
    output += "#SBATCH --ntasks=1"
    output += "\n"
    output += "#Set the time limit for the job - one hour is specified here"
    output += "\n"
    output += "#SBATCH --time=10:30:00"
    output += "\n"
    output += "#SBATCH --cpus-per-task=1"
    output += "\n"
    output += "# Load any environmental modules needed"
    output += "module load Python3"
    output += "module load gnu9"
    
    output += "# Move to the directory needed - defaults to the submission directory"
    output += "cd $HOME/research/Virgo/magphysParallel/output/"
    output += "\n"
    output += "# perform calculation"
    output += "\n"
    output += "\n"
    output += f"python {HOME}/github/virgosed/python/run1magphys.py {dirname}"

    outfname = f"JOB_{filename}.sh"
    outfile = open(outfname,'w')
    
    outfile.write(output)
    outfile.close()

    if submit:
        pass
        
