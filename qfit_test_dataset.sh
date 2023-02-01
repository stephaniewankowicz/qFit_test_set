y!/bin/bash
#$ -l h_vmem=4G
#$ -l mem_free=4G
#$ -t 1-174
#$ -l h_rt=2:00:00
#$ -R yes
#$ -V


#__________________SOURCE PHENIX/QFIT________________________________________________#
source /wynton/group/fraser/swankowicz/phenix-installer-1.19.2-4158-intel-linux-2.6-x86_64-centos6/phenix-1.19.2-4158/phenix_env.sh
export PATH="/wynton/home/fraserlab/swankowicz/anaconda3/bin:$PATH"
source activate qfit_debug

#________________PDB INFO__________________________________
PDB_file=/wynton/group/fraser/swankowicz/qfit_summit/expanded_dataset/list_of_pdbs_final.txt
PDB_dir='/wynton/group/fraser/swankowicz/qfit_summit/expanded_dataset/'
output_dir='/wynton/group/fraser/swankowicz/qfit_summit/expanded_dataset/output_data/single_structure/'


category='add_bulk'
cd ${output_dir}
PDB=$(cat $PDB_file | head -n $SGE_TASK_ID | tail -n 1)
echo ${PDB}


#______________________PDB STATS_______________________
#phenix.mtz.dump ${PDB_dir}/${PDB}/${PDB}.mtz > ${PDB_dir}/${PDB}/${PDB}.dump
#SPACE1=$(grep "^Space group symbol from file:" ${PDB_dir}/${PDB}/${PDB}.dump | awk '{print $6,$7}')
#UNIT1=$(grep "Unit cell:" ${PDB_dir}/${PDB}/${PDB}.dump | tail -n 1 | sed "s/[(),]//g" | awk '{print $3,$4,$5,$6,$7,$8}')
#RESO1=$(grep "^Resolution" ${PDB_dir}/${PDB}/${PDB}.dump | head -n 1 | awk '{print $4}')
#echo $PDB $RESO1 $SPACE1 $UNIT1 >> ${output_dir}/testset_space_unit_reso.txt

#find_largest_ligand.py ${PDB_dir}/${PDB}/${PDB}_qFit.pdb ${PDB}
b_factor.py ${PDB_dir}/${PDB}/${PDB}.pdb --pdb=${PDB}
phenix.rotalyze model=${PDB_dir}/${PDB}/${PDB}.pdb outliers_only=False > ${output_dir}/${PDB}_rotamer_output.txt

#b_fac=$(b_factor.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb --pdb=${PDB})
#phenix.rotalyze model=${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb outliers_only=False > ${output_dir}/${PDB}_rotamer_output.txt
#python /wynton/group/fraser/swankowicz/script/refine_log_parse.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.log ${PDB_dir}/${PDB}/${PDB}.updated_refine_001.log ${PDB}
