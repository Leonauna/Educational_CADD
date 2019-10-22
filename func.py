from rdkit.Chem import AllChem
import rdkit.Chem as Chem
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from rdkit.Chem import rdDepictor
import nglview as nv

def draw_pretty_pics(mol):
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

def show_3D(mol):
    """
    Given an RDKit mol object, return an NGLView 3D viewer
    """
    if mol.GetNumConformers() < 1:
        ids = AllChem.EmbedMultipleConfs(mol, numConfs=1)
    return nv.show_rdkit(mol)


