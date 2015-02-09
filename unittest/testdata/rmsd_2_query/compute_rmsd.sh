#!/bin/sh 

# bragg_ccp4_setup

from_pos="53"
to_pos="61"

while [ $to_pos -lt 121 ]
  do
    echo "FIT RESIDUE CA "$from_pos" TO "$to_pos > ./lsqkab.param
    echo "MATCH RESIDUE CA 2 TO 10" >> ./lsqkab.param
    new_pos=$[$from_pos-52]_1
    srmsd.py -p ./lsqkab.param -r ./3chy_a_2_1_query.pdb -d ./1ctf_a.pdb -l ./log.log 
    echo $new_pos
    # | awk -v value=$new_pos '{ print substr($1,1,6), new_pos, $2 }'
    from_pos=$[$from_pos+1]
    to_pos=$[$to_pos+1]
done

from_pos="2"
to_pos="10"

while [ $to_pos -lt 130 ]
  do
    echo "FIT RESIDUE CA "$from_pos" TO "$to_pos > ./lsqkab.param
    echo "MATCH RESIDUE CA 2 TO 10" >> ./lsqkab.param
    new_pos=$[$from_pos-1]_1
    srmsd.py -p ./lsqkab.param -r ./3chy_a_2_1_query.pdb -d ./3chy_a.pdb -l ./log.log 
    echo $new_pos
    # | awk -v value=$new_pos '{ print substr($1,1,6), new_pos, $2 }'
    from_pos=$[$from_pos+1]
    to_pos=$[$to_pos+1]
done

