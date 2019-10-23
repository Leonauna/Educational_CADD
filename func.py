from rdkit.Chem import AllChem
import rdkit.Chem as Chem
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from rdkit.Chem import rdDepictor
import nglview as nv
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from sklearn.metrics import r2_score
from rdkit import DataStructs
import multiprocessing

def view_estrogen_receptor():
    """
    load receptor structure and display with ngl
    """
    view = nv.show_pdbid("1g50", default=False)
    view._set_size("800px","400px")
    view.add_cartoon(":B or :C")
    view.add_ball_and_stick("(:B or :C) and EST")
    view.center(":B or :C")
    return view

def load_molecule_from_smiles(smiles):
    return Chem.MolFromSmiles(smiles)

def draw_pretty_2D_pics(mol):
    """
    Given an RDKit mol object, return an SVG picture
    """
    rdDepictor.Compute2DCoords(mol)
    mc_mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=True)
    drawer = rdMolDraw2D.MolDraw2DSVG(400, 200)
    drawer.DrawMolecule(mc_mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace("svg:", "")
    return SVG(svg)

def show_molecule_3D_structure(mol):
    """
    Given an RDKit mol object, return an NGLView 3D viewer
    """
    mol.RemoveAllConformers()
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=1)
    return nv.show_rdkit(mol)

def Xyfromdf(df, return_y):
    """Generate X (design matrix) and y (labels) for training a machine learning model from a dataframe of 
    structures. The structures are encoded using their Morgan fingerprint.

    Args:
        df: Pandas DataFrame with columns "CMPD_CHEMBLID", "STANDARD_VALUE", "CANONICAL_SMILES"
        return_y: whether or not to return labels y.

    Returns:
        2-tuple (X,y) where X is a dataframe and y is a series if return_y is True. Otherwise just X
    """
    # generate fingeprints: Morgan fingerprint with radius 2
	fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smile), 2) for smile in df["CANONICAL_SMILES"]]

	# convert the RDKit explicit vectors into numpy arrays
	np_fps = [np.zeros((1,)) for fp in fps]
	for i, fp in enumerate(fps):
	    DataStructs.ConvertToNumpyArray(fp, np_fps[i])
	X = pd.DataFrame(np.array(np_fps))
	if return_y:
		y = np.log(df["STANDARD_VALUE"])
		assert y.isna().sum()==0
		return X, y
	else:
		return X

def clean_data(df):
    """Remove missing values and duplicate CHEMBL IDs from a dataframe.

    Args:
        df: Pandas DataFrame with columns "CMPD_CHEMBLID", "STANDARD_VALUE", "CANONICAL_SMILES"

    Returns:
        The cleaned dataframe
    """
	#df = pd.read_csv("Data/training_data_raw.csv")
	df = df.dropna() # Remove missing values
	df = df.drop_duplicates(subset=["CMPD_CHEMBLID"])

	assert np.all([v==0 for k, v in {name : df[name].isna().sum() for name in df.columns}.items()])

	return df

def train_random_forest(
	training_data,
	seed=123,
	n_search_iter=10,
	k=4,
	n_jobs=multiprocessing.cpu_count()-1):
    """Train a random forest regression model on training data. The model predicts log(IC50) from smiley
    structure (encoded using Morgan fingerprints) 

    Args:
        training_data: a Pandas DataFrame with columns "CMPD_CHEMBLID", "STANDARD_VALUE", "CANONICAL_SMILES"
        seed: the seed for the pseudo-random number generator
        n_search_iter: the number of models in the random search cross-validation
        k: the number of folds in the cross-validation
        n_jobs: the number of cores used in the cross-validation

    Returns:
        The best performing random forest regression model according to the cross-validation.
    """
    np.random.seed(seed)
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
    """Predict the affinity for unseen structures using a trained random forest model.

    Args:
        random_forest: the random forest model
        test_data: a dataframe with columns "CMPD_CHEMBLID", "CANONICAL_SMILES"

    Returns:
        A dataframe with columns "CMPD_CHEMBLID", CANONICAL_SMILES", "predicted_affinity", where
        "predicted_affinity" is the random forest prediction of the raw IC50 value
    """
    
	X = Xyfromdf(test_data, False)
	test_data["predicted_affinity"] = np.exp(random_forest.predict(X))
	return test_data.sort_values(by="predicted_affinity", ascending=True)

