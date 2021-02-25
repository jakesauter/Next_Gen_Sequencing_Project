#! /bin/bash -l

#SBATCH  --partition=angsd_class
#SBATCH  --nodes=1
#SBATCH  --ntasks=1
#SBATCH  --job-name=project_data_download
#SBATCH  --time=12:00:00    # HH/MM/SS
#SBATCH  --mem=8G 
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

spack load sra-toolkit@2.10.7%gcc@6.3.0

data_dir="/athena/angsd/scratch/jns4001/project_data"

cd $data_dir

wget https://raw.githubusercontent.com/jakesauter/Next_Gen_Sequencing_Project/main/accession_tables/Sra_run_table.txt

cell_type="myeloid"
acc_nums=$(cat Sra_run_table.txt | grep $cell_type | cut -d',' -f1)

tempdir=/scratchLocal/${cell_type}_sra_download
mkdir $tempdir
cd $tempdir

for acc_num in $acc_nums ; do 
  fasterq-dump --outdir $cell_type $acc_num --temp=$tempdir
  gzip ${data_dir}/${cell_type}/${acc_num}.fastq
done
