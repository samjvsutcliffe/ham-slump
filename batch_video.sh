#!/bin/bash

# Request resources:
#SBATCH -c 8     # 1 entire node
#SBATCH --time=0:20:0  # 6 hours (hours:minutes:seconds)
#SBATCH --mem=4G      # 1 GB RAM
#SBATCH -p shared


module load python
module load ffmpeg
pip install pandas
pip install vtk

python make_video.py
