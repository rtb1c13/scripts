#!/bin/bash

# Script to reimage Tinker archive files into PDB for visualisation etc

# Options
prefix="19NT_dibut_run"          # Traj file prefix
center_atom=4              # 19NT carbonyl O
trajnum=2                  # No. of trajectories

# Cleanup
rm ${prefix}_all.arc xyz1.tmp

# Create xyzedit inputs
cat > xyz1.tmp << EOF
13
$center_atom

EOF
cat > xyz2.tmp << EOF
16

EOF

# Cat trajectories
for ((i=1;i<=$trajnum;i++))
do
 gunzip ${prefix}${i}.arc.gz
 cat ${prefix}${i}.arc >> ${prefix}_all.arc
 gzip ${prefix}${i}.arc
done

# Center, reimage and change to pdb (Make sure .seq file is present if needed!)
xyzedit ${prefix}_all.arc -k ${prefix}.key < xyz1.tmp
xyzedit ${prefix}_all.xyz -k ${prefix}.key < xyz2.tmp && mv ${prefix}_all.xyz_2 ${prefix}_all_reimaged.arc

if [ -f ${prefix}1.seq ]
then
 echo "Seq file present - using this to define topology"
 cp ${prefix}1.seq ${prefix}_all_reimaged.seq
else
 echo "Seq file NOT present - topology will be undefined!"
fi

xyzpdb ${prefix}_all_reimaged.arc -k ${prefix}.key
sed -i '/CONECT/d' ${prefix}_all_reimaged.pdb

# Final cleanup
rm ${prefix}_all.xyz ${prefix}_all.arc
gzip ${prefix}_all_reimaged.arc
echo "Done"
