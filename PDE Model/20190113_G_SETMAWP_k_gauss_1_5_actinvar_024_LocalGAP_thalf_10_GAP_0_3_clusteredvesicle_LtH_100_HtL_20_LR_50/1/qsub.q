#!/bin/tcsh
#
#SBATCH -o screen.log
#SBATCH -e screen.err
#SBATCH  --mem=1G

/opt/apps/MATLAB/R2012b/bin/matlab -nosplash -singleCompThread -r 'RUN;quit'
