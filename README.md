# Model Training

directory: `train_TestSet`

Create a Conda environment NN-gpu for training the model using the file ` NeuralNet.yml `

Create a folder `CellType` for the cell type. Move the  `.h5`  training and validation files to `CellType/data` as `dataset-train.h5` and `dataset-val.h5`respectively. The dataset files for training, validation and testing can be generated using the Nextflow Pipeline below.

For example the files `data.type66.rep2-train.rnd.h5` and  `data.type66.rep2-val.rnd.h5` are moved to `type66/data/{dataset-train.h5,dataset-val.h5}`.

Run the following command `scriptGPU.sh train type66 3`, preferably with GPU access.


# Nextflow Pipeline

Using the DeepTact pipeline we want to predict the relation between a promoter-enhancer pair using input features corresponding to Chromatin Openness and training labels generated from promoter capture Hi-C data. The input features for the Chromatin Openness score can be generated from either ATAC-seq/Dnase-seq data. 

The pipeline is built such that it can be modified easily to create multiple output datasets with lego-like modularity for input datasets. To that end we have seperated the pipeline into 5 steps. Steps 1 and 2 generate only species specific features and can be completely generated given the reference genome. Step 3 and Step 4 require .bam and ChiCAGO processed .tsv files as input, respectively. Step 5 combines all the steps above to generate hdf5 dataset for input to the Neural Network.

Run the pipeline inside an environment with nextflow, & python. To run the test script above, run the nextflow pipeilne in the 5 `nxf_TestSet/5-TestSet_nxf/run_*.sh` scripts. The input dataset files are the following:
- Dataset paths for reference genome files in `pr-enh-prep.config` in `params.genomes`
- Dataset paths for the aligned .bam ATAC-seq/DNase-seq files in `co-score-prep.config` in `params.bamInput`
- Dataset paths for the ChICAGO processes PCHiC input files in `pchic-prep.config` in `params.hic_input`

Modify the following .config files to change the output dataset parameters : `nxf_TestSet/{1,2,3,4,5}-TestSet_nxf/conf/*.config`. To modify the nextflow computation options modify the parameters in : `nxf_TestSet/{1,2,3,4,5}-TestSet_nxf/nextflow.config`.


## Step1: Promoter and Enhancer Lists

Pipeline produces the species specific .bed files with lists of both regulatory elements, with lengths `promoter_window/enhancer_window + augment_length` & `bgWindow` : promoters and enhancers. It also generates a .bed file containing comprehensive list of all promoter-enhancer pairs within `max_dist_pr_enh` distance in all chrommosomes.


directory: `nxf_TestSet/1-TestSet_nxf`  
config file(s): `conf/pr-enh-prep.config`  
run: ` ./run_pr_enh.sh -resume `  

### Parameters

filename : `pr-enh-prep.config`

Sample p **dataDir**: path to species specific datasets -- used inside `params.genomes`
* **refgen** : reference genome for processing  : mm10, hg18, hg38
* **bgWindow**: length of background window for normalization and in Chromatin Openness Score calculation (default: 1e6)
* **augment_length**: length of augmentation of promoter and enhancer window for data augmentation (default: 1000)
* **augment_step**: step size for moving window for data augmentation (default: 50)
* **promoter_window, enhancer_window**: promoter/enhancer window size (default: 1000 & 2000 respectively)
* **max_dist_pr_enh**: max distance between brute force generated promoter-enhancer pairs  (default: 1e6)
* **all_chrom**: chromosome list for the selected species
* **species**['\<key\>']: : chromosome list for the species corresponding to *\<key\>*  
* **genomes**: Dataset files for reference genomes  
► *key*: reference genome name  
► *fasta*        : fasta file  
► *fasta_idx*    : fasta index file  
► *promoter_gtf* : .gtf file with a list of promoters  
► *enhancer_bed* : .bed file with a list of promoters  
► *enhBed_headers* : headers in the enhancer .bed file  

Example:
```  genomes {
    'mm10' {
      fasta        = "$dataDir/mm10/mm10.fa"
      fasta_idx    = "$dataDir/mm10/mm10.fa.fai"
      promoter_gtf = "$dataDir/mm10/gencode/mm10_gencodevM25.gtf.gz"
      enhancer_bed = "$dataDir/enhancers_fantom5/F5.mm10.enhancers.bed"
      enhBed_headers = "['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']"
    }
  }

  species {
    'hg' {
      all_chrom = "['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']"
    }
    'mm' {
      all_chrom = "['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY']"
    }
  }

  all_chrom = params.species[params.refgen.substring(0,2)].all_chrom
```



## Step2:  Genomic Sequence Feature

Pipeline produces the species specific .npz files containing the one hot encoded genomic sequence for the promoter and enhancer elements generated in step 1.


directory: `nxf_TestSet/2-TestSet_nxf`  
config file(s): `conf/genome-seq-prep.config`, `conf/pr-enh-prep.config`  
run: ` ./run_genome_seq.sh -resume `  

### Parameters: 

filename : `genome-seq-prep.config`

* **store_dir**: base directory with processed input files  required from output of step-1
* **regelem_inp**: path to file specifying the output from step 1.  

Header format for **regelem_inp**:  
    - `regelem` : regultory element name (enhancer/promoter/brutePairs)  
    - `main`: bed file path relative to `params.store_dir`  
    - `background`: background bed  path relative to `params.store_dir`  


## Step3: Chromatin Openness Feature

Pipeline produces the gzipped hdf5 files (.h5.gz) for the Chromatin Openness Score (CO score) feature from the ATAC-seq data for the promoter and enhancer lists from Step 1.

directory: `nxf_TestSet/3-TestSet_nxf`  
config file(s): `conf/co-score-prep.config`, `conf/pr-enh-prep.config`  
run: ` ./run_co_score.sh --dev -resume `  

### Parameters: 

filename : `co-score-prep.config`

* **bamInput**: paths for the .bam files for the ATAC-seq data
* **store_dir**: base directory with processed input files  required from output of step-1 (similar to Step 2)
* **regelem_inp**: path to file specifying the output from step 1 (similar to Step 2)


Header format for **bamInput** :  
`celltype`: cell type name  
`repetition`: repetition number  
`bam`: absolute path for .bam file  

## Step4: Positive Interaction Labels

Pipeline produces the .csv files from PCHiC data to use as input for Neural Network training/testing. If the options are set for cross-validation splitting, the dataset is split into test, val, training, and left-overs datasets by chromosome taking care to match as closely the ratio specified in `params.split`.

directory: `nxf_TestSet/4-TestSet_nxf`  
config file(s): `conf/pchic-prep.config`, `conf/pr-enh-prep.config`  
run: ` ./run_pchic.sh --dev -resume `  

### Parameters: 

filename : `co-score-prep.config`

* **store_dir**: base directory with processed input files  required from output of Step 1, 2 and 3 (similar to Step 2)
* **hic_input**: "proc_files/hic_input.csv"
* **cellTypes**: list for cell types in the ChiCAGO processes PCHiC data
* **dnaseq_inp**: path to file specifying the output from step 2
* **regelem_inp**: path to file specifying the output from step 1 (similar to Step 2)
* PCHi-C data processing options
    * **pos_threshold**: ChiCAGO score threshold for a PCHiC interaction to be considered positive. (Threshold for negative interactions is taken to be 0)
    * **hic_split_process**: split hic tsv into `hic_split_process` parts for further parallel processing after separation by chromosome 
    * NeuralNet Training Data input params:
        * **hic_augment_factor**: data augmentation factor for PCHiC data
        * **neg_pos_interac_ratio**: `'inf' / positive integer (1,2,3,..)`
        * **split**: data split ratio for cross-validation

Parameters to generate the features for all promoter-enhancer pairs without data augmentation:
```
  hic_augment_factor = 1
  neg_pos_interac_ratio = 'inf'
  split = []//all
```

Sample parameters to generate the features for all promoter-enhancer pairs with data augmentation:
```      
    hic_augment_factor = 20 
    neg_pos_interac_ratio = 4 
    split = [0.2, 0.2] // test / val / train
```

Header format for **hic_input**:  
    - `hic_input_name`: dataset name  
    - `path`: asbolute path to processed ChiCAGO output in .tsv format  

Header format for **dnaseq_inp**:  
    - `enhancer`: path relative to `params.store_dir` to enhancer npz file from step 2  
    - `promoter`: path relative to `params.store_dir` to promoter npz file from step 2  


## Step5: hdf5 Dataset for NeuralNet

Pipeline produces the hdf5 files (.h5) containg the features and training data to use as input for Neural Network training/testing corresponding to the input PCHiC file (from output of step 4) . The output hdf5 includes additional info such as regulatory element location, name.

directory: `nxf_TestSet/5-TestSet_nxf`  
config file(s): `conf/feature-prep.config`, `conf/pr-enh-prep.config`  
run: ` ./run_feature_gen.sh --dev -resume  --randomize`  

### Parameters: 

filename : `feature-prep.config`

* **regelem_inp**: path to file specifying the output from step 1 (similar to Step 2)
* **coscore_inp**: path to file specifying the output from step 3
* **dnaseq_inp**: path to file specifying the output from step 2 (similar to Step 4)
* **chunks_val, chunks_test, chunks_train, chunks_all**: chunks for processing the validation, test, train, and all/comprehensive pr-enh pair data respectively
* **pchic_inp**: path to file specifying the output from step 4
* **randomize**: randomly sorts the output .hdf5 file when true -- for easy loading for NeuralNet training in `DeepTact.py`

Header format for **coscore_inp**:  
    - `celltype`: celltype name  
    - `repetition`: repetition number  
    - `regelem`: regulatory element name : enhancer/promoter  
    - `h5`: absolute path for .h5.gz file from step 3  

Header format for **pchic_inp**:  
    - `celltype`: celltype name. should be same as in `coscore_inp`  
    - `label`: the last characters of the label should match : `test` / `train` / `val` / `all` / `left-overs`. This is to use the `chunks_val`, `chunks_test`, `chunks_train`, `chunks_all` variable properly  
    - `hic`: absolute path for .csv file from step 4  





