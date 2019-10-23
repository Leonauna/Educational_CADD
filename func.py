from rdkit.Chem import AllChem
import rdkit.Chem as Chem
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from rdkit.Chem import rdDepictor
import nglview as nv

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
    if mol.GetNumConformers() < 1:
        ids = AllChem.EmbedMultipleConfs(mol, numConfs=1)
    return nv.show_rdkit(mol)


