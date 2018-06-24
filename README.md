bacpacs version 0.0.5<br>User's guide
==============================

Overview
--------

bacpacs is a bacterial pathogenicity classification python module, based on the following paper: "BacPaCS – Bacterial Pathogenicity Classification via Sparse-SVM", by Eran Barash, Neta Sal-Man, Sivan Sabato, and Michal Ziv-Ukelson (Submitted). It can be used for classification using a pre-trained model, or for generating a new model from labeled training data of sequenced proteomes. The training pipeline:

1. bacpacs selects the 10% longest protein out of the set of all proteins from all training samples.

2. Using CD-HIT, bacpacs clusters the selected proteins. This results in
   clusters, or protein families (PFs). Each PF is represented by its
   longest protein.

3. For each organism in the training set, bacpacs compares the
   organism's proteins to the representatives of the PFs, using CD-HIT-2D.
   This results in a binary feature vector for each organism. The
   vector indices represent the different PFs.
   If a cell with index _i_ has a value `True`, the organism has a protein
   similar to the representative of PF _i_. Otherwise, the cell's value is `False`.

4. The binary feature vectors of the organisms in the training set (and their known
   pathogenicity labels), can then be used to train a linear SVM model (using l1
   norm as penalty). Other models can be used as well. 


Installation and dependencies
-----------------------------

Bacpacs can be installed via one of the following two alternatives:
1. Run in a linux terminal: `$ pip install --index-url https://test.pypi.org/simple/ bacpacs` (recommended)
2. Clone or download bacpacs Github [repository](https://github.com/barashe/bacpacs.git) and run `pip install -e path/to/bacpacs`

Dependencies:

-   Operating system: Linux. CD-HIT only runs on Linux OS, and so does bacpacs.
-   CD-HIT: bacpacs requires CD-HIT to be installed. It is also recommended to add CD-HIT to the PATH variable. CD-HIT can be downloaded from its official [website](http://weizhongli-lab.org/cd-hit/).

-   Python 2.7.

-   Several python packages, including numpy, scipy, and scikit-learn. If
    you are missing a package, you will get an error stating which package
    you are missing. Then you can simply install the package using
    [pip](https://pypi.python.org/pypi/pip/) or
    [EasyInstall](https://wiki.python.org/moin/EasyInstall).


Running
-------
Below are elaborated running examples of the two possible bacpacs schemes:
1. Predicting data using bacpacs pre-trained model: bacpacs comes with a pre-trained model, used in the bacpacs paper. The pre-trained model can be easily downloaded and used.
2. Training and using a model: bacpacs can also be used to translate a training set of organisms into a feature vector that can be fed into an SVM training module, and to then translate a test set into a feature vector which uses the same features as the training set. The pathogenicity of the organisms in the test set can then be predicted using the trained model. The organisms in both the training set and the test set are fed as raw amino acid fasta files (faa files).

The example code below should be used in Python 2.7. Full documentation of each of the methods appears in the code. This example code can be found in [examples](https://github.com/barashe/bacpacs/tree/master/examples).



## Predicting data using bacpacs pre-trained model

### Imports

In[1]
```python
import bacpacs
from sklearn.svm import LinearSVC
```

The bacpacs pre-trained model can be found in https://github.com/barashe/bacpacs/tree/master/trained. It can also be downloaded and loaded directly from Python as described below.

### Load the pre-trained model

In[2]
```python
bp, clf = bacpacs.load_trained_model(output_dir='out1')
```
Out[2]

    Retrieving files from github
    Downloading full_bacpacs.pkl
    Downloading linearsvc_full.pkl
    Downloading protein_families


load_trained_model() returns a Bacpacs object, and a sklearn.svm.LinearSVC trained classifier. We will need bp for feature extraction, and clf for the prediction itself. 
Here 'bacpacs' is an empty folder, and 'out1' is created when invoking load_trained_model(). These names can be set according to preference. An existing output_dir is acceptable as well, but beware of file overrun. 

### Download toy data
Use your real data, or download bacpacs toy data (also available on [Github](https://github.com/barashe/bacpacs/blob/master/toy.tar.gz)):

In[3]
```python
bacpacs.download_toy_data('.')
```
Out[3]

    Toy data stored in ./toy

You can move genomes from toy/train to toy/validate and vice versa, or even delete some. The destination folder can be set by replacing '.' with the desired destination.

### Extract features for test organisms
This step takes a while, depending on your machine:

In[4]
```python
bp.genomes_vs_pfs('toy/validate', n_jobs=0)
```
Out[4]

    Running genomes against protein families representatives
    Genome cluster files are stored in out1/pred_clusters


genomes_vs_pfs() uses CD-HIT-2D to run all test genomes against the pre-established protein families. Assigning n_jobs=0 tells CD-HIT-2D to use all available CPUs for each genome. However, the genomes are processed sequentially independently. Running genomes_vs_pfs() on several machines can save plenty of time. 
<br><br>Next, extract the features and store them in variables:

In[5]

```python
X_pred = bp.extract_features(feats_type='pred')
```

### Read pathogenicity labels from a csv file
bacpacs.read_labels() takes a csv file with two columns and no headers: the first column should list genome ids, and the second column should list pathogenicity labels. Genome ids should match the original genome file names. For example: for a genome file named org1.faa, the csv file should list an 'org1' genome id. Pathogenicity should be boolean: True for pathogens, False for non-pathogens. Note that we include the corresponding feature matrix (X_pred), to ensure that the order of the returned labels corresponds to the order of the genomes in the feature matrix.

In[6]
```python
y_true = bacpacs.read_labels('toy/labels.csv', X=X_pred)
```
### Predict pathogenicity

In[7]
```python
y_pred = clf.predict(X_pred)
```
Compute accuracy:

In[8]
```python
(y_pred == y_true).mean()
```

Out[8]


    0.80000000000000004



Or simply:

In[9]

```python
clf.score(X_pred, y_true)
```
Out[9]



    0.80000000000000004



## Training and using a model

### Imports

In[1]
```python
import bacpacs
from sklearn.svm import LinearSVC
```

### Download toy data
Use your real data, or download the bacpacs toy data (also available on [Github](https://github.com/barashe/bacpacs/blob/master/toy.tar.gz)):


In[2]
```python
bacpacs.download_toy_data('.')
```
Out[2]

    Toy data stored in ./toy


### Initialize a Bacpacs object

In[3]
```python
bp = bacpacs.Bacpacs('out2')
```

### Merge training data

Merge all training .faa files into on large .faa file.

In[4]
```python
bp.merge_genome_files('toy/train/', output_path=None)
```
Specify an output path to override the default destination directory.

Out[4]

    Saving merged proteins as out/merged.faa


Note that it is possible to specify a different output path.

### Reduce the merged file to the 10% longest proteins

In[5]
```python
bp.reduce(long_percent=10, merged_path=None, output_path=None)
```

Out[5]

    Saving reduced proteins as out/reduced.faa

long_percent can be set to any value between 1 and 100. The number of selected proteins is rounded down if the requested percentage does not result in a whole number. 
The default merged_path is <output_directory>/merged.faa and the default output_path is <output_directory>/reduced.faa. Both can be set using the appropriate arguments of 'reduce'.


### Cluster the training proteins into protein families

In[6]

```python
bp.create_pfs(memory=800, n_jobs=0, cdhit_path=None, reduced_path=None, output_path=None)
```

Out[6]

    Clustering genomes.
    cd-hit -i out/reduced.faa -o out/protein_families -c 0.4 -n 2 -M 800 -T 0
    Clustering finished successfully. Protein families dumped in out/protein_families


If CD-HIT is included in the system's path, there is no need to provide a path to 'cdhit_path'. If cd-hit is not included in the path, a valid path must be provided. 
Note that we are using n_jobs=0, to use all available CPUs. 

### Extract features

In[7]
```python
bp.genomes_vs_pfs('toy/train/', feats_type='train', n_jobs=0)
```

Out[7]

    Running genomes against protein families representatives
    Genome cluster files are stored in out/train_clusters

In[8]

```python
bp.genomes_vs_pfs('toy/validate/', feats_type='pred', n_jobs=0)
```

Out[8]

    Running genomes against protein families representatives
    Genome cluster files are stored in out/pred_clusters


In[9]

```python
X_train = bp.extract_features(feats_type='train')
X_pred = bp.extract_features(feats_type='pred')
```

### Read pathogenicity labels from a csv file
bacpacs.read_labels() takes a csv file with two columns and no headers: the first column should list genome ids, and the second column should list pathogenicity labels. Genome ids should match the original genome file names. For example: for a genome file named org1.faa, the csv file should list an 'org1' genome id. Pathogenicity should be boolean: True for pathogens, False for non-pathogens. Note that we include the corresponding feature matrices (X_train, X_pred), to ensure that the orders of the returned label sets correspond to the orders of the genomes in the feature matrices.

In[10]

```python
y_train = bacpacs.read_labels('toy/labels.csv', X_train)
y_pred = bacpacs.read_labels('toy/labels.csv', X_pred)
```

Load a linear SVM object:

In[11]

```python
clf = LinearSVC(penalty='l1', dual=False)
```

Fit it with the created features and labels:

In[12]
```python
clf.fit(X_train, y_train)
```


Out[12]

    LinearSVC(C=1.0, class_weight=None, dual=False, fit_intercept=True,
         intercept_scaling=1, loss='squared_hinge', max_iter=1000,
         multi_class='ovr', penalty='l1', random_state=None, tol=0.0001,
         verbose=0)



Compute the accuracy pf the model's prediction:

In[13]

```python
clf.score(X_pred, y_pred)
```


Out[13]

    0.80000000000000004


