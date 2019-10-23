import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from sklearn.metrics import r2_score
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import multiprocessing

def Xyfromdf(df, return_y):
    # generate fingeprints: Morgan fingerprint with radius 2
	fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smile), 2) for smile in df["CANONICAL_SMILES"]]

	# convert the RDKit explicit vectors into numpy arrays
	np_fps = [np.zeros((1,)) for fp in fps]
	for i, fp in enumerate(fps):
	    DataStructs.ConvertToNumpyArray(fp, np_fps[i])
	X = pd.DataFrame(np.array(np_fps))
	if return_y:
		y = pd.Series(np.log(df["STANDARD_VALUE"].values))
		assert y.isna().sum()==0
		return X, y
	else:
		return X

def clean_data(df):
	#df = pd.read_csv("Data/training_data_raw.csv")
	df = df.dropna() # Remove missing values
	df = df.drop_duplicates(subset=["CMPD_CHEMBLID"])

	assert np.all([v==0 for k, v in {name : df[name].isna().sum() for name in df.columns}.items()])

	return df

def train_random_forest(
	training_data,
	test_size=0.3,
	seed=123,
	n_search_iter=10,
	k=4,
	n_jobs=multiprocessing.cpu_count()-1):

	X, y = Xyfromdf(training_data, True)

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
	).fit(X, y)

	return rf.best_estimator_

def predict_affinity(
	random_forest,
	test_data):
    
	X = Xyfromdf(test_data, False)
	test_data["predicted_affinity"] = random_forest.predict(X)
	return test_data.sort_values(by="predicted_affinity", ascending=True)