# Model Training
dir: ```train_TestSet```

Create a Conda environment NN-gpu for training the model using the file ``` NeuralNet.yml ```

Create a folder ```CellType``` for the cell type. Move the  ```.h5```  training and validation files to ```CellType/data``` as ```dataset-train.h5``` and ```dataset-val.h5```respectively. The dataset files for training, validation and testing can be generated using the Nextflow Pipeline below.

For example the files ```data.type66.rep2-train.rnd.h5``` and  ```data.type66.rep2-val.rnd.h5``` are moved to ```type66/data/{dataset-train.h5,dataset-val.h5}```.

Run the following command ```scriptGPU.sh train type66 3```, preferably with GPU access.


# Nextflow Pipeline

dir: ```nxf_TestSet```

## Step1: 

dir: ```1-TestSet_nxf```

config file(s): ```conf/pr-enh-prep.config```

run: ``` ./run_pr_enh.sh -resume ```


## Step2: 

dir: ```2-TestSet_nxf```

config file(s): 
```conf/genome-seq-prep.config```, ```conf/pr-enh-prep.config```

run: ``` ./run_genome_seq.sh -resume ```

## Step3: 

dir: ```3-TestSet_nxf```

config file(s): ```conf/co-score-prep.config```, ```conf/pr-enh-prep.config```

run: ``` ./run_co_score.sh --dev -resume ```

## Step4: 

dir: ```4-TestSet_nxf```

config file(s): ```conf/pchic-prep.config```, ```conf/pr-enh-prep.config```

run: ``` ./run_pchic.sh --dev -resume ```

## Step5: 

dir: ```5-TestSet_nxf```

config file(s): ```conf/feature-prep.config```, ```conf/pr-enh-prep.config```

run: ``` ./run_feature_gen.sh -resume  --randomize```
