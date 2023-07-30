# RNNHC
A hybrid error correction algorithm based on Recurrent Neural Network
## Environment requirements
Make sure the server is installed Minimap2, Samtools, and Anaconda.  
The environment requirement: python3.6, Tensorflow 2.3  

## Correction steps
Run prepare correction_data_pb.sh to correct a pacbio long reads, the same execution as ONT reads  
the script includes the process of evaluation.  

## Overview of the files
mp_SR_as_label.py is used to encode the pileup file  
tf_train_for_SR_as_label.py is used to train the model  
tf_infer_for_SR_as_label.py is used to predict and decode the corrected reads  
coverage_count2.0.py is used to count the coverage of the long reads
## 
for yeast PacBio dataset, minimap2 should add the parameters "minimap2 --split-prefix=tmp -ax sr -t"
