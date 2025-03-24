#!/usr/bin/env bash
#SBATCH --job-name=install-seurat
#SBATCH --output=/cfs/earth/scratch/paleslui/BATH/output/install_seurat_%j.log
#SBATCH --error=/cfs/earth/scratch/paleslui/BATH/output/install_seurat_%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=earth-4
#SBATCH --constraint=rhel8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --mail-user=paleslui@students.zhaw.ch
#SBATCH --mail-type=ALL

module purge
module load USS/2022
module load gcc/9.4.0-pe5.34

export PATH=/cfs/software/uss/2022/spack/linux-rocky8-x86_64/gcc-9.4.0/miniconda3-4.12.0-7rii.../bin:$PATH

# Remove old environment if exists
conda env remove --name seurat_env -y

# Create Seurat environment
conda env create -f /cfs/earth/scratch/paleslui/BATH/SLURMs/enviroment.seurat.yml

echo "Seurat environment setup completed!"
