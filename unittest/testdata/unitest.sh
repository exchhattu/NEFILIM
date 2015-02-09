#!/bin/sh

# echo "Make backup"
# if [ -f ./1ctf_a_32_1.out ]; then 
#   mv ./1ctf_a_32_1.out ./1ctf_a_32_1_b.out
# fi
# 
# if [ -f ./1ctf_a_40_1.out ]; then 
#   mv ./1ctf_a_40_1.out ./1ctf_a_40_1_b.out
# fi
# 
# if [ -f ./3chy_a_15_1.out ]; then 
#   mv ./3chy_a_15_1.out ./3chy_a_15_1_b.out
# fi

../../TestCheck -r @./flags

# echo "Removed output files."
# rm -rf ./1ctf_a_32_1.out 
# rm -rf ./1ctf_a_40_1.out 
# rm -rf ./3chy_a_15_1.out 
# 
# echo "Backup"
# mv ./1ctf_a_32_1_b.out ./1ctf_a_32_1.out  
# mv ./1ctf_a_40_1_b.out ./1ctf_a_40_1.out 
# mv ./3chy_a_15_1_b.out ./3chy_a_15_1.out 

