#!/bin/bash
# $1 is the first command line argument which is the number of replicas.
if [ $# -ne 1 ];then
    echo "Incorrect no of arguments"
    exit
else
    num_reps=$1
fi

  dir_data=$PWD
  mkdir $dir_data/1-traj_analysis
  dir_out=$dir_data/1-traj_analysis
  for ((n=1; n<=num_reps; n++)); do
    dir_in=$dir_data/rep$n
    awk '($1=="BOND"||$1=="VDWAALS"||$1=="1-4")' $dir_in/_MMPBSA_complex_pb.mdout.all | cut -c11-24,36-49,64-77 | awk -f $SCRIPTS_PATH/t-ser-mmpb.awk > $dir_out/mmpb_com.dat
    awk '($1=="BOND"||$1=="VDWAALS"||$1=="1-4")' $dir_in/_MMPBSA_complex_gb.mdout.all | cut -c11-24,36-49,64-77 | awk -f $SCRIPTS_PATH/t-ser-mmpb.awk > $dir_out/mmgb_com.dat
    awk '($1 != "#Frame"){print $2*0.0072}' $dir_in/_MMPBSA_complex_gb_surf.dat.all > $dir_out/esurf_com.dat
    awk '($1=="BOND"||$1=="VDWAALS"||$1=="1-4")' $dir_in/_MMPBSA_receptor_pb.mdout.all | cut -c11-24,36-49,64-77 | awk -f $SCRIPTS_PATH/t-ser-mmpb.awk > $dir_out/mmpb_rec.dat
    awk '($1=="BOND"||$1=="VDWAALS"||$1=="1-4")' $dir_in/_MMPBSA_receptor_gb.mdout.all | cut -c11-24,36-49,64-77 | awk -f $SCRIPTS_PATH/t-ser-mmpb.awk > $dir_out/mmgb_rec.dat
    awk '($1 != "#Frame"){print $2*0.0072}' $dir_in/_MMPBSA_receptor_gb_surf.dat.all > $dir_out/esurf_rec.dat
    awk '($1=="BOND"||$1=="VDWAALS"||$1=="1-4")' $dir_in/_MMPBSA_ligand_pb.mdout.all | cut -c11-24,36-49,64-77 | awk -f $SCRIPTS_PATH/t-ser-mmpb.awk > $dir_out/mmpb_lig.dat
    awk '($1=="BOND"||$1=="VDWAALS"||$1=="1-4")' $dir_in/_MMPBSA_ligand_gb.mdout.all | cut -c11-24,36-49,64-77 | awk -f $SCRIPTS_PATH/t-ser-mmpb.awk > $dir_out/mmgb_lig.dat
    awk '($1 != "#Frame"){print $2*0.0072}' $dir_in/_MMPBSA_ligand_gb_surf.dat.all > $dir_out/esurf_lig.dat

    paste $dir_out/mmpb_com.dat $dir_out/esurf_com.dat | awk -f $SCRIPTS_PATH/t-ser-mmpbsa.awk > $dir_out/rep$n-pb-com.dat 
    paste $dir_out/mmpb_rec.dat $dir_out/esurf_rec.dat | awk -f $SCRIPTS_PATH/t-ser-mmpbsa.awk > $dir_out/rep$n-pb-rec.dat 
    paste $dir_out/mmpb_lig.dat $dir_out/esurf_lig.dat | awk -f $SCRIPTS_PATH/t-ser-mmpbsa.awk > $dir_out/rep$n-pb-lig.dat 
    paste $dir_out/mmgb_com.dat $dir_out/esurf_com.dat | awk -f $SCRIPTS_PATH/t-ser-mmgbsa.awk > $dir_out/rep$n-gb-com.dat 
    paste $dir_out/mmgb_rec.dat $dir_out/esurf_rec.dat | awk -f $SCRIPTS_PATH/t-ser-mmgbsa.awk > $dir_out/rep$n-gb-rec.dat 
    paste $dir_out/mmgb_lig.dat $dir_out/esurf_lig.dat | awk -f $SCRIPTS_PATH/t-ser-mmgbsa.awk > $dir_out/rep$n-gb-lig.dat 
  done

rm $dir_out/mm?b_???.dat $dir_out/esurf_???.dat
