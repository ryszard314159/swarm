#!/bin/bash
#-----------------------------------------------------------------------------
# rocs_score1.sh: Generates a ROCS score report file, given a template mol file 
# and a file of SMILES strings corresponding to molecules to be overlaid
#
# Input: requires 2 commandline args:
#   - arg1 must be the name of the SMILES filename of the mol to be scored
#   - arg2 must be the filename of the template mol - the ligand/substrate to be matched to (with 3D coords)
# Output: a ROCS report file with the same name as the input SMILES file, but '.rpt' extension
#
# This script is also used as the ROCS driver by SWARM
#-----------------------------------------------------------------------------

# Clip suffix from SMILES filename and generate other filenames from resulting basename
baseName=${1%.smi}
flipperOut=$baseName.ism
omegaOut=$baseName.sdf
rocsOut=$baseName.rpt

# Get process ID and use that as the base name for a temporary file for the dumb Openeye logo output
processID=$$
oeLogoOut=$$.tmp

source /usr/contrib/etc/bash.modules
module load python
module load oeprod

# Try sending the OpenEye terminal trace to a file so we can watch SWARM trace better:
flipper -in $1 -out $flipperOut >& $oeLogoOut
omega2 -in $flipperOut -out $omegaOut  >& $oeLogoOut
rocs -query $2 -dbase $omegaOut -rankby combo -nostructs true -besthits 0 -reportfile $rocsOut >& $oeLogoOut

# cleanup intermediate files; keep log files for debugging failures:
rm -f $flipperOut
rm -f $omegaOut
rm -f $oeLogoOut

#rm -f omega2.log
rm -f omega2.parm
rm -f omega2.pending.ism
rm -f omega2_status.txt
#rm -f rocs.log
rm -f rocs.parm
rm -f rocs_1.status
