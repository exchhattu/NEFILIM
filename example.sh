#!/bin/sh

#at source the rosetta library since it is not in ~/.bashrc

echo "give the rosetta database path and redirect it to the flags" 
echo "$ROSETTA_PATH/rosetta_database" > flags 

./NEFILIM -r @flags \
         --path_to_scores ./inputs/apmds.lst \
         --path_to_models ./inputs/model.lst \
         --path_to_native ./inputs/1ctf_atom.pdb \ # just for computing rmsd to fragments
         --window_size 9 \
         --cutoff 0.30 \
         --proportion 1 \
         --number_of_fragment 25 \
         --fragment_per_cluster 5 \
         --number_of_template 5 \
         --cluster_radius 1.00

mv ./N_AP_RDS_fragment.fragK ./N_AP_RDS_fragment.frag9
mv ./N_AP_RDS_template.fragK ./N_AP_RDS_template.frag9
mv ./N_SLIDING_MEAN_AP_RDS_fragment.fragK ./N_SLIDING_MEAN_AP_RDS_fragment.frag9
mv ./N_SLIDING_MEAN_AP_RDS_template.fragK ./N_SLIDING_MEAN_AP_RDS_template.frag9
mv ./rsd_scores.out ./rsd_scores_frag9.out

echo "for three residues "
./NEFILIM -r @flags \
         --path_to_scores ./inputs/apmds.lst \
         --path_to_models ./inputs/model.lst \
         --path_to_native ./inputs/1ctf_atom.pdb \
         --window_size 3 \
         --cutoff 0.30 \
         --proportion 1 \
         --number_of_fragment 25 \
         --fragment_per_cluster 5 \
         --number_of_template 5 \
         --cluster_radius 0.20

mv ./N_AP_RDS_fragment.fragK ./N_AP_RDS_fragment.frag3
mv ./N_AP_RDS_template.fragK ./N_AP_RDS_template.frag3
mv ./N_SLIDING_MEAN_AP_RDS_fragment.fragK ./N_SLIDING_MEAN_AP_RDS_fragment.frag3
mv ./N_SLIDING_MEAN_AP_RDS_template.fragK ./N_SLIDING_MEAN_AP_RDS_template.frag3
mv ./rsd_scores.out ./rsd_scores_frag3.out
