## Table of Contents
1. [Run on gitpod](#run-on-gitpod)
2. [Run on HPC](#run-on-hpc)
3. [Run Locally](#run-locally)
4. [Pipeline Details](#pipeline-details)
    - [Profiles](#profiles)
    - [Configuration Files](#configuration-files)
    - [SLURM Configuration](#slurm-configuration)

## Run on gitpod
This method involves uploading the required input files directly to gitpod and running the nextflow pipeline from there. It is the fastest and most efficient way to execute the pipeline.

#### **1. Open gitpod**
Click the link below to launch gitpod with the repository:  

[Open in gitpod](https://gitpod.io/new/#https://github.com/acg-team/Bioinfo4B)  

#### **2. Navigate to the pipeline folder**
Once gitpod loads, navigate to the following directory:  

```
cd Genomics/nextflow_bioinformatics4beginners
```

This is where the nextflow pipeline is located.  

#### **3. Download input files**  
We need to create an `input` folder and place the required files there. These files are not hosted on GitHub due to their size. You can download them from [google drive](https://drive.google.com/file/d/13_op58XI3L3S2DGsd1r7zLmh9GYyevXF/view?usp=sharing).

After you upload them on gitpod, you should have following structure
```

 Genomics
 │ ├── nextflow_bioinformatics4beginners
 │ │ ├── input
 │ │ │   ├── reads.fastq
 │ │ │   ├── chrM.fa
 │ │ ├── main.nf
 │ │ ├── (other pipeline files)

```

#### **4. Run the pipeline**
After placing the input files, execute the Nextflow pipeline with the following command:  

```bash
nextflow run main.nf -profile docker \
    --input $(realpath input/reads.fastq)  \
    --reference $(realpath input/chrM.fa) \
    --outdir $(realpath results)
```

This will process the provided sequencing data and store the results in the `results` folder.  

## Run on HPC

#### **1. Set up working directory**
Navigate to your scratch space and create a folder for the pipeline:
```sh
cd $LSFM_CLUSTER_SCRATCH_USER_PATH
mkdir bio4beginners_nextflow
cd bio4beginners_nextflow
```
#### **2. Load required modules**
Load the necessary HPC modules:
```sh
module load USS/2022 gcc/9.4.0-pe5.34 miniconda3/4.12.0
```
If you are using conda for the first time, you need to initialize it:
```sh
conda init bash
```
After running this command, **close and reopen your terminal**. After that you will need to start from step 1.

#### **3. Create conda environment**
First, set up bioconda according to the bioconda documentation, notably setting up channels:
```sh
conda config --add channels bioconda
conda config --add channels conda-forge
```

Now, create the conda environment:
```sh
conda create --name env_nf nextflow
```

And activate it:
```sh
conda activate env_nf
```

#### **4. Copy source code**
Copy the pipeline source folder and prepare the output directory:
```sh
cp -r /cfs/earth/scratch/shared/bioinfo4beginners/Genomics/bio4beginners_nextflow/source .
cd source
mkdir results
```

#### **5. Run the pipeline**
Execute the Nextflow pipeline using test data:
```sh
nextflow run main.nf -profile conda -c slurm.config \
    --input $(realpath input/reads.fastq)  \
    --reference $(realpath input/chrM.fa) \
    --outdir $(realpath results)
```
This command:
- Uses the **conda profile** to manage dependencies.
- Runs with **SLURM** using `slurm.config`.
- Converts all file paths to **absolute paths** to avoid path issues.

After execution, results will be stored in the `results/` directory.

## Run locally

#### **1. Download source code and input files**

Download `.zip` file from [google drive](https://drive.google.com/file/d/1A3vcVaQmviO27aJFSsh2r4mkBbFXPKTY/view?usp=sharing) and navigate to the project directory.

#### **2. Create conda environment**

First, set up bioconda according to the bioconda documentation:

```sh
conda config --add channels bioconda
conda config --add channels conda-forge
```

Now, create the conda environment:

```sh
conda create --name env_nf nextflow
```

and activate it:

```sh
conda activate env_nf
```

### **3. Running the pipeline locally** 

To execute the pipeline, navigate to the **source directory** and run:

```sh
nextflow run main.nf -profile conda \
    --input $(realpath input/reads.fastq)  \
    --reference $(realpath input/chrM.fa) \
    --outdir $(realpath results)
```

### **Execution steps**
1. **conda environment setup**  
   - The pipeline will first create a **new conda environment** from `Env_Genomics.yml`.  
   - This step may take approximately **10 minutes**.  

2. **Pipeline execution**  
   - After the environment is ready, the execution should take about **5 minutes**.  
   - Results will be stored in the `results/` directory.  


### **Running with a container**  
If conda fails, you can run the pipeline using **containers** (e.g., docker, podman, or singularity).  

For example, for docker container, switch to the **docker profile** by running:

```sh
nextflow run main.nf -profile docker \
    --input $(realpath input/reads.fastq)  \
    --reference $(realpath input/chrM.fa) \
    --outdir $(realpath results)
```

If your container runtime requires **root privileges**, you must execute nextflow with `sudo`:

```sh
sudo nextflow run main.nf -profile docker \
    --input $(realpath input/reads.fastq)  \
    --reference $(realpath input/chrM.fa) \
    --outdir $(realpath results)
```


## Pipeline details

![Pipeline Diagram](img/mermaid-diagram-2025-03-10-194316.png)

### **Profiles**  

Profiles define the **execution environment** for Nextflow and are specified using `-profile <profile_name>`.  This pipeline supports `conda`, `docker`, `podman`, `singularity`

#### **conda profile**  
The **conda environment** is specified inside `nextflow.config` using `process.conda`:  

1. **Default**  
In `nextflow.config`, `process.conda` is set to:  
```
process.conda = "$projectDir/Env_Genomics.yml"
```
This creates a new conda environment from `Env_Genomics.yml` when the pipeline runs.  

2. **HPC execution**  
In `slurm.config`, `process.conda` is set to:  
```nextflow
process.conda = "/cfs/earth/scratch/shared/bioinfo4beginners/Genomics/Env_Genomics"
```
On Earth cluster, an existing conda environment is already available, so the pipeline directly uses it without creating a new environment.  

### Configuration files

Each Nextflow pipeline has a default configuration file called `nextflow.config`. This file allows defining parameters to be used in `main.nf`.

All input parameters (e.g., input files, reference genome, output directory) can be specified in `nextflow.config` using `params`. These parameters are accessible in `main.nf` through `params.<param_name>`.

For example, instead of specifying inputs in the command line, they can be provided in `nextflow.config`:

```nextflow
params {
    input = "input/reads.fastq"
    reference = "input/chrM.fa"
    outdir = "results"
}
```

This allows running the pipeline without explicitly passing parameters in the command line:

```sh
nextflow run main.nf -profile conda
```

If additional configuration files are needed (e.g., for running the pipeline on an HPC system), they can be included using the `-c` parameter when executing Nextflow.

For example:

```sh
nextflow run main.nf -profile conda -c slurm.config
```

This allows different configurations to be applied without modifying the default `nextflow.config` file.


### SLURM configuration
The `slurm.config` file is used to run Nextflow with SLURM. It:
- Submits each process as a SLURM job
- Defines resource limits (CPU, memory, partition, etc.)
- Allows customization per process

Typically, to run a script on an HPC system with SLURM, we use `sbatch`. However, Nextflow supports SLURM as an executor and can submit jobs automatically. To enable this, set `process.executor` to `slurm` in the configuration file.

For Earth HPC, it is also required to specify the partition and resource limits. These details can be configured in `slurm.config`. It is also possible to set specific resource requirements for individual processes within this file.
