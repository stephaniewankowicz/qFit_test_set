#!/bin/bash
#$ -l h_vmem=24G
#$ -l mem_free=24G
#$ -t 1-174
#$ -l h_rt=30:00:00
#$ -pe smp 4
#$ -R yes
#$ -V
#this script will run qfit based on the input PDB names you have.

#________________________________________________INPUTS________________________________________________#
#base_dir='/wynton/group/fraser/swankowicz/200104_macrodomain/'  #where the folders are located
#PDB_file=/wynton/group/fraser/swankowicz/200104_macrodomain/refine_again_221107.txt #list of PDB IDs
#PDB_file=/wynton/group/fraser/swankowicz/PDZ/pdb_list.txt
#base_dir='/wynton/group/fraser/swankowicz/PDZ/'

PDB_file=/wynton/group/fraser/swankowicz/qfit_summit/expanded_dataset/list_of_pdbs_final.txt
base_dir='/wynton/group/fraser/swankowicz/qfit_summit/expanded_dataset/'
category='add_bulk'

#PDB_file=/wynton/group/fraser/swankowicz/qfit_summit/expanded_cryoEM/pdb.txt
#base_dir='/wynton/group/fraser/swankowicz/qfit_summit/expanded_cryoEM/'
#PDB_file=/wynton/group/fraser/swankowicz/designed_proteins/pdbs.txt
#base_dir='/wynton/group/fraser/swankowicz/designed_proteins/'

#PDB_file=/wynton/group/fraser/swankowicz/COVID19/ACE2/pdb.txt
#base_dir='/wynton/group/fraser/swankowicz/COVID19/ACE2/'

#base_dir='/wynton/group/fraser/swankowicz/qfit_debug/fall_summit_2022/qfit_stability/b_sampling/' #seg_4' #base folder (where you want to put folders/pdb files
#PDB_file=/wynton/group/fraser/swankowicz/qfit_debug/fall_summit_2022/qfit_stability/pdbs.txt 

export OMP_NUM_THREADS=1

#________________________________________________SET PATHS________________________________________________#
source /wynton/group/fraser/swankowicz/phenix-installer-1.19.2-4158-intel-linux-2.6-x86_64-centos6/phenix-1.19.2-4158/phenix_env.sh
export PATH="/wynton/home/fraserlab/swankowicz/anaconda3/bin:$PATH"
source activate qfit_debug
which python

#________________________________________________RUN QFIT________________________________________________#
 PDB_dir=$(cat $PDB_file | head -n $SGE_TASK_ID | tail -n 1)
 #echo $PDB_dir 
 #if [[ -z "$TMPDIR" ]]; then
 #  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
 #  mkdir -p "$TMPDIR"
 #  export TMPDIR
 #fi

 #cd ${TMPDIR}

 #cp -R ${base_dir}/${PDB_dir}/ ${TMPDIR}/
 cd ${base_dir}/${PDB_dir}
 PDB=$(echo ${PDB_dir:0:4})
 echo ${PDB}
 mkdir ${category}
 cd ${category}
 pwd

#__________________________________DETERMINE FOBS v IOBS v FP__________________________________
mtzmetadata=`phenix.mtz.dump "${PDB}.mtz"`

# List of Fo types we will check for
obstypes="FP FOBS F-obs I IOBS I-obs FC"

# Get amplitude fields
ampfields=`grep -E "amplitude|intensity" <<< "${mtzmetadata}"`
ampfields=`echo "${ampfields}" | awk '{$1=$1};1' | cut -d " " -f 1`

# Clear xray_data_labels variable
xray_data_labels=""

for field in ${ampfields}; do
  echo $field
  # Check field in obstypes
  if grep -F -q -w $field <<< "${obstypes}"; then
    echo "found obs"
    # Check SIGFo is in the mtz too!
    if grep -F -q -w "SIG$field" <<< "${mtzmetadata}"; then
      xray_data_labels="${field},SIG${field}";
      break
    fi
  fi
done


#if [[ -e segment_multiconformer_model2.pdb ]]; then 
#   echo 'done'
#else
   #phenix.composite_omit_map ${PDB}.mtz ${PDB}.pdb omit-type=refine nproc=16 r_free_flags.generate=True
#   qfit_protein composite_omit_map.mtz -l 2FOFCWT,PH2FOFCWT ${PDB}.updated_refine_001.pdb -p 8 -d ${base_dir}/${PDB_dir}/${category}
#fi

if [[ -e ${PDB}_qFit.pdb ]]; then 
   echo 'done'
else
   cp ../segment_constraints/qFit_occupancy.params . 
   cp ../multiconformer_model2.pdb . 
   cp ../${PDB}.mtz .
   qfit_final_refine_xray.sh ${PDB}.mtz multiconformer_model2.pdb
fi
#qfit_final_refine_xray.sh ${PDB}_updated.pdb.updated_refine_001.mtz multiconformer_model2.pdb


#cp -R ${TMPDIR}/${PDB_dir}/ ${base_dir}
#done
#done < $PDB_file
