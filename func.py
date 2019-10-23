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
	test_data["predicted_affinity"] = np.exp(random_forest.predict(X))
	return test_data.sort_values(by="predicted_affinity", ascending=True)

