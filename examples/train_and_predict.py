# Training and using a model

# Imports
import bacpacs
from sklearn.svm import LinearSVC

# Download toy data
bacpacs.download_toy_data('.')

# Initialize a Bacpacs object
bp = bacpacs.Bacpacs('out')

# Merge training data
bp.merge_genome_files('toy/train/', output_path=None)

# Reduce the merged file to the 10% longest proteins
bp.reduce(long_percent=10, merged_path=None, output_path=None)

# Cluster training proteins to protein families
bp.create_pfs(memory=800, n_jobs=0, cdhit_path=None, reduced_path=None, output_path=None)

# Extract features
bp.genomes_vs_pfs('toy/train/', feats_type='train', n_jobs=0)
bp.genomes_vs_pfs('toy/validate/', feats_type='pred', n_jobs=0)

X_train = bp.extract_features(feats_type='train')
X_pred = bp.extract_features(feats_type='pred')

# Read pathogenicity labels from a cv file
y_train = bacpacs.read_labels('toy/labels.csv', X_train)
y_pred = bacpacs.read_labels('toy/labels.csv', X_pred)

# Initialize and fit a classifier
clf = LinearSVC(penalty='l1', dual=False)
clf.fit(X_train, y_train)

# Score your prediction
clf.score(X_pred, y_pred)

