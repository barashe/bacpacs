bacpacs V 0.0.5 -- User's guide
==============================

Overview
--------

bacpacs is a bacterial pathogenicity classification python module. bacpacs
trains using sequenced proteomes. The training pipeline:

1. Considering all proteins, from all training samples, bacpacs samples
   the 10% longest proteins.

2. Using CD-HIT, bacpacs clusters the sampled proteins. This results in
   clusters, or protein families (PFs). Each PF is represented by its
   longest protein.

3. For each organism in the training set, bacpacs compares the
   organism's proteins to the clusters' representatives, using CD-HIT-2D.
   This results in a binary feature vector for each organism. The
   vectors indices, represent the different clusters' representatives.
   If a cell with index _i_ has a value `True`, the organism has a protein
   similar to representative _i_. Otherwise, the cell's value is `False`.

4. The training set feature vectors (and their known
   pathogenicity labels), can then be used to train a Linear SVM model (using l1
   norm as penalty). Other models can be used as well. 


Installation and dependencies
-----------------------------

Two alternatives to install bacpacs:
1. pip install --index-url https://test.pypi.org/simple/ bacpacs (recommended)
2. clone or download bacpacs Github [repository] (https://github.com/barashe/bacpacs.git) and run `pip install -e path/to/bacpacs`

Dependencies:

-   CD-HIT -- bacpacs requires CD-HIT to be installed. It is also recomended to add CD-HIT to the PATH variable. CD-HIT can be downloaded from its official [website](http://weizhongli-lab.org/cd-hit/).
-   CD-HIT only runs on Linux OS, and so does bacpacs. 

-   Python 2.7.

-   Several python packages, including numpy, scipy, scikit-learn. If
    you are missing a package, you'll get an error stating which package
    you are missing. Then you can simply install the package using
    [pip](https://pypi.python.org/pypi/pip/) or
    [EasyInstall](https://wiki.python.org/moin/EasyInstall).


Running
-------
Below are elaborated running examples of the two possible bacpcas schemes:
1. Predicting data using bacpacs pre-trained model - bacpacs comes with a pre-trained model, used in the bacpacs paper (in publication process). The module enables simple download and usage of the trained model.
2. Training and using a model - bacpacs can also be used to create features to train a new model, from raw amino acid fasta files (faa files).

### Imports


```python
import bacpacs
from sklearn.svm import LinearSVC
```

## Predicting data using bacpacs pre-trained model

The bacpacs trained model can be found in https://github.com/barashe/bacpacs/tree/master/trained. However, a simpler way to download and load the model is described here.

### Load the pre-trained model


```python
bp, clf = bacpacs.load_trained_model(output_dir='out1')
```

    Retrieving files from github
    Downloading full_bacpacs.pkl
    Downloading linearsvc_full.pkl
    Downloading protein_families


load_trained_model() returns a Bacpacs object, and a sklearn.svm.LinearSVC trained classifier. We will need bp for feature extraction, and clf for the prediction itself. 
Here 'bacpacs' is an empty folder, and 'out1' is created when invoking load_trained_model(). An existing output_dir is acceptible as well, but beware of file overrun. 

### Download toy data
Use your real data, or download bacpacs toy data (also availbale on [Github](https://github.com/barashe/bacpacs/blob/master/toy.tar.gz)):

```python
bacpacs.download_toy_data('.')
```

    Toy data stored in ./toy

You can move genomes from toy/train to toy/validate and vice versa, or even delete some.

### Extract prediction features
This step takes a while, depending on your machine:


```python
bp.genomes_vs_pfs('toy/validate', n_jobs=0)
```

    Running genomes against protein families representatives
    Genome cluster files are stored in out1/pred_clusters


genomes_vs_pfs() uses CD-HIT-2D to run all prediction genomes against the pre-established protein families. Assigning n_jobs=0, tells CD-HIT-2D to use all available CPUs. However, each genome is processed at a time, independently. Running extract_features() on several machines, can save plenty of time. 
<br><br>Next, we'll extract the features and store them in variables:


```python
X_pred = bp.extract_features(feats_type='pred')
```

### Read pathogenicity labels from a csv file
bacpacs.read_labels() takes a csv file, which is expected to have to columns (no headers): genome id, and pathogenicity. Genome ids should match the original genomes file names. For example: for a genome file named org1.faa, the csv file should include an 'org1' genome id. Pathogenicity should be boolean: True for pathogens, False for non-pathogens. Note that we include the corresponding feature matrix (X_pred), to insure the correct genomes and their order.

```python
y_true = bacpacs.read_labels('toy/labels.csv', X=X_pred)
```

X is an optional parameter which uses X_pred to make sure X_pred and y_pred are in the same order.<br><br>
We can now use clf for prediction:


```python
y_pred = clf.predict(X_pred)
```
Compute accuracy:
```python
(y_pred == y_true).mean()
```




    0.80000000000000004



Or simply:


```python
clf.score(X_pred, y_true)
```




    0.80000000000000004



## Training and using a model

### Download toy data
Use your real data, or download bacpacs toy data (also availbale on [Github](https://github.com/barashe/bacpacs/blob/master/toy.tar.gz)):

```python
bacpacs.download_toy_data('.')
```

    Toy data stored in ./toy


### Initialize a Bacpacs object


```python
bp = bacpacs.Bacpacs('out2')
```

### Merge training data

We begin by merging all training .faa files to one large .faa file.


```python
bp.merge_genome_files('toy/train/', output_path=None)
```

    Saving merged proteins as out/merged.faa


Note that it is possible to specify a different output path.

### Reduce the merged file to the 10% longest proteins


```python
bp.reduce(long_percent=10, merged_path=None, output_path=None)
```

    Saving reduced proteins as out/reduced.faa


### Cluster training proteins to protein families


```python
bp.create_pfs(memory=800, n_jobs=0, cdhit_path=None, reduced_path=None, output_path=None)
```

    Clustering genomes.
    cd-hit -i out/reduced.faa -o out/protein_families -c 0.4 -n 2 -M 800 -T 0
    Clustering finished successfully. Protein families dumped in out/protein_families


Since cd-hit is included in path, there is no need to supply a path to 'cdhit_path'. If cd-hit is not included in the path, a valid path must be supplied. 
Note we are using n_jobs=0, to use all available CPUs. 

### Extract features


```python
bp.genomes_vs_pfs('toy/train/', feats_type='train', n_jobs=0)
```

    Running genomes against protein families representatives
    Genome cluster files are stored in out/train_clusters



```python
bp.genomes_vs_pfs('toy/validate/', feats_type='pred', n_jobs=0)
```

    Running genomes against protein families representatives
    Genome cluster files are stored in out/pred_clusters



```python
X_train = bp.extract_features(feats_type='train')
X_pred = bp.extract_features(feats_type='pred')
```

### Extract pathogenicity labels

### Read pathogenicity labels from a csv file
bacpacs.read_labels() takes a csv file, which is expected to have to columns (no headers): genome id, and pathogenicity. Genome ids should match the original genomes file names. For example: for a genome file named org1.faa, the csv file should include an 'org1' genome id. Pathogenicity should be boolean: True for pathogens, False for non-pathogens. Note that we include the corresponding feature matrices (X_train and X_pred), to insure the correct genomes and their order. This way, the same labels.csv file can be used for both sets.


```python
y_train = bacpacs.read_labels('toy/labels.csv', X_train)
y_pred = bacpacs.read_labels('toy/labels.csv', X_pred)
```

Load a linear SVM object:
```python
clf = LinearSVC(penalty='l1', dual=False)
```

Fit it with the fatures and labels created:
```python
clf.fit(X_train, y_train)
```




    LinearSVC(C=1.0, class_weight=None, dual=False, fit_intercept=True,
         intercept_scaling=1, loss='squared_hinge', max_iter=1000,
         multi_class='ovr', penalty='l1', random_state=None, tol=0.0001,
         verbose=0)



Compute prediction's accuracy:
```python
clf.score(X_pred, y_pred)
```




    0.80000000000000004


