### HPC Workload Management via SLURM

#### SLURM Script: Header and Directives
```
#!/bin/bash   # Use the Bash shell to run this script

#SBATCH --mail-user=user@mail.fhsu.edu  # Email address for job updates
#SBATCH --job-name=spades_assembly      # Descriptive job name
#SBATCH --mail-type=all                 # Send email on start, end, and failure
#SBATCH --mem=2048                      # Allocate 2048 MB of memory
#SBATCH --cpus-per-task=6               # Request 6 compute nodes
#SBATCH --time=12:00:00                 # Set maximum runtime to 12 hours
```

#### SLURM Script: Software Setup and Command Execution
```
# Load software module (ensures environment is set up correctly)
module load SPAdes

# Create an output directory if it does not exist
mkdir -p spades_assembly_out

# Run SPAdes genome assembler with paired-end reads
spades.py --isolate \
  -1 plastomeOnlyReads_R1.fastq \   # Forward reads file
  -2 plastomeOnlyReads_R2.fastq \   # Reverse reads file
  -o spades_assembly_out \          # Output directory
  -t 6                              # Use 6 threads (matches: #SBATCH --cpus-per-task=6)
```

#### Complete SLURM Script Example
```
#!/bin/bash
#SBATCH --mail-user=user@mail.fhsu.edu
#SBATCH --job-name=spades_assembly
#SBATCH --mail-type=all
#SBATCH --mem=2048
#SBATCH --cpus-per-task=6
#SBATCH --time=12:00:00

module load SPAdes            # Load software
mkdir -p spades_assembly_out  # Generate output directory

# Run spades
spades.py --isolate -1 plastomeOnlyReads_R1.fastq \
-2 plastomeOnlyReads_R2.fastq -o spades_assembly_out -t 6
```

#### Submitting and Monitoring Runs via SLURM
```
sbatch spades_plastome.sbatch  # submits the job

squeue -u $USER                # shows pending/running 
                               # jobs of user; shows job IDs
                         
scontrol show job 3708305      # displays detailed information 
                               # about a specific job (resources, 
                               # node allocation, time limits, etc.)

cat slurm-3708305.out          # SLURM output provides regular 
                               # standard output/error of job
```

#### Cancelling a Run Prematurely
```
scancel 3708305                # cancels particular job 
                               # (using job ID)
```
