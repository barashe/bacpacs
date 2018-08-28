# Predict data using bacpacs pre-trained model

# Imports
import bacpacs

# Load the pre-trained model
bp, clf = bacpacs.load_trained_model(output_dir='out1')

# Download toy data
bacpacs.download_toy_data('.')

# Extract features for test organisms
bp.genomes_vs_pfs('toy/validate', n_jobs=0, memory=102400)
X_pred = bp.extract_features(clusters_dir='out1/pred_clusters/', feats_type='pred')


# Read pathogenicity labels from a csv file
y_true = bacpacs.read_labels('toy/labels.csv', X=X_pred)

# Prediction
y_pred = clf.predict(X_pred)
print 'Prediction score: ', (y_pred == y_true).mean()

# Or simply:
print 'prediction score (quick method): ', clf.score(X_pred, y_true)
