import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from sklearn.metrics import r2_score
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import multiprocessing

def get_best_candidates(
	n_candidates=10,
	test_size=0.3,
	seed=123,
	n_search_iter=10,
	k=4,
	n_jobs=multiprocessing.cpu_count()-1):
	"""Returns a dataframe of CHEMBL IDs and SMILES structures of the most promising
	compounds for targeting estrogen receptor alpha according to a random forest
	regression model.

	Args:
		n_candidates; the number of candidates (default is 10)
		test_size: the fraction of the compounds used in the test set (default is 0.3)
		seed: the seed for the pseudo-random number generator (default is 123)
		n_search_iter: the number of models tested in the random search cross-validation (default is 10)
		k: the number of folds in the stratified cross-validation (default is 4)
		n_jobs: the number of threads used in the cross-validation (deafult is all the available cores minus one)

	Returns:
		A dataframe with two columns (CHEMBL ID, SMILEY structure) and n_candidate rows, containing the most promising
		compounds. The compounds are sorted such that the most and least promoising compounds are in the first and last
		rows respectively.
	"""
	print("Getting best candidates...", end="")
	np.random.seed(seed)

	df = pd.read_csv("Data/bioactivity_clean.csv")

	# generate fingeprints: Morgan fingerprint with radius 2
	fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smile), 2) for smile in df["CANONICAL_SMILES"]]

	# convert the RDKit explicit vectors into numpy arrays
	np_fps = [np.zeros((1,)) for fp in fps]
	for i, fp in enumerate(fps):
	    DataStructs.ConvertToNumpyArray(fp, np_fps[i])
	X = pd.DataFrame(np.array(np_fps))
	y = pd.Series(np.log(df["STANDARD_VALUE"].values))
	assert y.isna().sum()==0

	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=seed)

	#sns.distplot(y, rug=True)

	# Do cross-validation for random forest

	rf = RandomizedSearchCV(
	    RandomForestRegressor(),
	    {
	        'n_estimators' : np.arange(100, 1000, 100),
	        'max_features' : ['sqrt', 'log2'],
	        'max_depth' : [None] + list(range(100))
	    },
	    n_iter=n_search_iter,
	    cv=k,
	    random_state=seed,
	    n_jobs=n_jobs,
	    verbose=0
	).fit(X_train, y_train)

	# plt.scatter(y_test, rf.predict(X_test))
	# plt.xlabel("truth")
	# plt.ylabel("predicted")
	# plt.title("r^2 on test set: {:.2f}".format(rf.score(X_test, y_test)))

	# Get the best candidates
	y_test_hat = pd.DataFrame(rf.predict(X_test), index=X_test.index, columns=["yhat"])
	y_test_hat = y_test_hat.sort_values("yhat")

	if y_test.shape[0] < n_candidates:
		n_candidates = y_test.shape[0]

	best_candidates = df.iloc[ y_test_hat.index[:n_candidates] ][ ["CMPD_CHEMBLID", "CANONICAL_SMILES"] ] # These are the best candidates
	print("done!")
	return best_candidates