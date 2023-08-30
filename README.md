# RNNHC
A hybrid error correction algorithm for long reads based on Recurrent Neural Network
## Environment requirements
Make sure the server is installed Minimap2, Samtools, and Anaconda.  
The environment requirement: python3.6, Tensorflow-gpu 2.3  

## Correction steps
Run prepare correct_pac.sh to correct a pacbio long reads, the same execution as ONT reads  
Make sure the .pileup file is generated after "minimap2 alignment" and "samtools mpileup"  
the script includes the process of evaluation.  

## Overview of the files
mp_SR_as_label.py is used to encode the pileup file for train 
tf_train_for_SR_as_label.py is used to train the model  
tf_infer_for_SR_as_label.py is used to predict and decode the corrected reads  
coverage_count3.0.py is used to count the percent of bases in long reads that covered by short reads
## 
for yeast PacBio dataset, minimap2 should add the parameters "minimap2 --split-prefix=tmp -ax sr -t"
