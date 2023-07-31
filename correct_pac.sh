
# This .sh include aligned only the short reads to the long reads to be corrected.
# then merge the two mpileup, divied the mpileup by linecounts to adapt the multi preprocessing.
# generate a whole hdf for train and a set of split_hdf for infer
# usage: sh sr_as_label.sh long short ref


long_reads = $1
short_reads = $2
reference = $3

mkdir split_mpileup
mkdir split_hdf
mkdir split_fasta


# align and generate a complete mpileup

minimap2 --split-prefix=tmp -ax sr -t 30 $1 $2 -a --secondary=no -o short2long.bam 
samtools sort short2long.bam -o short2long.sorted.bam
samtools index short2long.sorted.bam
samtools mpileup short2long.sorted.bam -a -s -f $1 -o mpileup_genome.pileup -a

rm -rf *.bam *.fai *.bai


# split the mpileup into a folder
python ./split_file.py --mpileup mpileup_genome.pileup --output-folder ./split_mpileup
rm -rf ./split_mpileup/manifest

# encode the splited mpileup into a complete HDF file and write the uncovered sequences out uncovered.fasta for train
python ./mp_SR_as_label.py --mpileup-folder ./split_mpileup --output-file covered_fasta.hdf --uncovered-fasta uncover.fasta

# encode the splited mpileup into splited HDF files for infer
python ./split_mpileup_for_SR_as_label.py --mpileup-folder ./split_mpileup --output-folder ./split_hdf

# train
python ./tf_train_for_SR_as_label.py --hdf-files covered_fasta.hdf --check-point-path ./weight/cp.ckpt --model_hdf5_path ./weight/model.h5

# infer
python ./tf_infer_for_SR_as_label.py --hdf-folder ./split_hdf --output-folder ./split_fasta --model-hdf5-path ./weight/model.h5


cat split_fasta/*fasta > corrected.fasta
cat uncover.fasta corrected.fasta > rnnhc_"$1"

# evaluate the results

sh ./evaluate.sh $3 rnnhc_"$1" pb

# sh sr_as_label.sh F.antasticus_long_error.fa F.antasticus_short.fa F.antasticus_genome.fa


rm -rf split_mpileup
rm -rf split_hdf






