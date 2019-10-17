from rdkit import Chem
from rdkit.Chem import MCS
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdchem as rdchem


mol1 = Chem.MolFromSmiles("Cc1ccccc1")
mol2 = Chem.MolFromSmiles( "CCc1ccccc1" )
mol3 = Chem.MolFromSmiles( "Oc1ccccc1" )
mol4 = Chem.MolFromSmiles( "COc1ccccc1" )



import re
import numpy as np # argsort


def read_xyz():
    return(11)

def read_sdf():
    return(22)

def read_mol():
    return(33)

def read_no():
    return(44)

def read_file(file):
    suffix = re.search('\.[^\.]*$',file).group(0)
    return {
        '.xyz': read_xyz(file),
        '.mol': read_mol(file),
        '.sdf': read_sdf(file),
    }.get(suffix, read_no)


# Fracture query_mol into fragments consisting of core and sidechains
def FragmentMol(query_mol, query_mcs_indices):
  answer = rdchem.EditableMol(query_mol)
  all_query_bonds = query_mol.GetBonds()
  for bond_index in range(len(query_mol.GetBonds()):
    begin_atom_idx = all_query_bonds[bond_index].GetBeginAtom().GetIdx()
    end_atom_idx = all_query_bonds[bond_index].GetEndAtom().GetIdx()
    # if the bond is between mcs core and a side-chain - remove it
    rdchem.EditableMol.RemoveBond(answer, begin_atom_idx, end_atom_idx) if ( bool(begin_atom_idx in query_mcs_indices) != bool(end_atom_idx in query_mcs_indices) ) else None
  return(answer.GetMol())
    

def MapIndexQ2T(query_index, query_mcs_indices, template_mcs_indices):
  # get index of query_index in query_mcs_indices
  # return template_mcs_indices[index]
  # else return -1


class Sidechain:
  def __init__():
    anchor = (-1, -1, -1) # path of length 3 starting in atom where sidechain connects that goes through 2 other atoms that are part of a skeleton
      # these are stored not as indices in context of query mol but instead as order in query_mcs_indices
      # it is assumed this order is identical in template_mcs_indices and thus is used for accessing atom indices directly in template_mol
    members = list() # contains indices of atoms (in context of whole molecule) that belong to this sidechain
    build_order = list() # indices of members (in context of sidechain) listed in ascending order of distance from skeleton contact; ties are undefined but irrelevan
    zpath = list() # list of triplets - first three indices (in context of whole molecule) of path to skeleton contact
    zvals = list() # list of triplets - distance/angle/dihedral to atoms from zpath
    dists_to_anchor = list() # distance of shortest path to skeleton contact - used for filling build_order
    atom_type = list()
  
  # add atom to sidechain by adding its index to members, first three steps of a path to anchor for z-matrix like reconstruction,
  # zvals for distance+angle+dihedral and total distance to anchor
  def AddAtom(query_mol, index):
  
  # fill anchor values by providing index of a skeleton contact
  def Anchor(query_mol, skeleton_contact, query_mcs_indices):
    for mcs_idx_it in query_mcs_indices:
      path = Chem.rdmolops.GetShortestPath(query_mol, skeleton_contact, mcs_idx_it)
      if ( len(path) == 3 & path[1] in query_mcs_indices & path[2] in query_mcs_indices ):
        anchor[0] = np.where(np.array(path[0]) == query_mcs_indices) # skeleton_contact
        anchor[1] = np.where(np.array(path[1]) == query_mcs_indices) # way_to_anchor
        anchor[2] = np.where(np.array(path[2]) == query_mcs_indices) # anchor
        return(0)
    return(1)
  
  # fill build_order
  def OrderUp(query_mol):
    for it in members:
      path = Chem.rdmolops.GetShortestPath(query_mol, it, self.skeleton_contact) # path will include starting and ending index
      dists_to_anchor.add(len(path)) # remember the length
      # we pay special attention if the path is not long enough to build the atom from zvals of this sidechain only
      if (len(path)) == 2 ):
        # atom is adjacent to skeleton contact
        # dihedral and angle from template is used; only distance is taken from query
      if (len(path)) == 3 ):
        # atom is 2 steps away from skeleton contact
        # only dihedral from template is used - IS THIS WELL DEFINED? LAST ATOM OF DIHEDRAL IS FROM QUERY, IS DIFF FROM SKELETON CONTACT AND IS AMBIGUOUS. DOES THIS MATTER?
        # !! this should be ok if angles are not changed
        # angle and  distance are taken from query
      if (len(path)) > 3 ):
        # all zvals take from sidechain
        zpath.add(path[-3:]) # add last three steps of the path
        # measure distance it - path[1] and use it as zval[0]
        # measure angle it - path[1] - path[2] and use it as zval[1]
        # measure dihedral it - path[1] - path[2] - path[3] and use it as zval[2]
    # now time to order
    build_order = np.argsort(dists_to_anchor)

# end of Sidechain class

class MorphMol:
  # rdkit molecule
  # sidechain list

  def Reconstruct():
  # reconstruct atom Idx using coordinates of 3 zpath reference atoms 
  # These are rdkit.Geom.Point3D ??? - pick something with already implemented measurements of distance/angles/dihedrals



  
    




# xyz to Mol
## call babel
## babel -ixyz

query_mol = readfile(file1)
template_mol = readfile(file2)

# mcs
skeleton_mcs = MCS.FindMCS([query_mol, template_mol])
skeleton_mol = Chem.MolFromSmarts(skeleton_mcs.smarts)


min_overlap = 6  # Require the overlap to be at least min_overlap atoms

if ( len(skeleton_mol.GetAtoms()) < min_overlap):
    print("These molecules share less than min_overlap atoms. Quitting.")


query_mcs_matches = query_mol.GetSubstructMatches(skeleton_mol)
print("query_mol contains overlapping fragment " + len(query_mcs_matches) + " times.")
template_mcs_matches = template_mol.GetSubstructMatches(skeleton_mol)
print("template_mol contains overlapping fragment " + len(template_mcs_matches) + " times.")

for query_mcs_match_index in range(len(query_mcs_matches)):
  for template_mcs_match_index in range(len(template_mcs_matches)):
    # get indices of mcs match
    query_mcs_indices = query_mcs_matches[query_mcs_match_index]
    query_mcs_indices_set = set(query_mcs_indices)
    template_mcs_indices = template_mcs_matches[template_mcs_match_index]

    # create a copy of mol and remove all bonds between mcs and side_chains
    fragmented_query_mol = FragmentMol(query_mol, query_mcs_indices)

    # Get sets of indices of individual sidechains  
    fragments_tuples = Chem.rdmolops.GetMolFrags(fragmented_query_mol)
    fragments_sets = {}
    for i in fragments_tuples: fragments_sets.add(set(i))
    sidechains_sets = {}
    for i in fragments_sets:
      if ( i != query_mcs_indices_set ): sidechains_sets.add(i)

    # are there any cycles that include both a sidechain and more than one skeleton atom? If so, such morph is invalid
    # we check this by checking neighbours of all atoms in a sidechain and adding those that are part of skeleton to a set. If len(set) > 1 => reject morph

    all_sidechains = list()

    for sidechain_it in sidechains_sets:
      skeleton_contact_current = set()
      for atom_it in sidechain_it:
        for neighbours_it in query_mol.GetAtoms()[atom_it].GetNeighbours():
          atom_idx = neighbours_it.GetIdx()
          skeleton_contact_current.add(atom_idx) if atom_idx in query_mcs_indices else None
      
      if ( len(skeleton_contact) == 0 ): # throw a problem
      if ( len(skeleton_contact) > 1 ): # throw a problem
      if ( len(skeleton_contact) == 1 ):
        new_sidechain = SideChain()
        new_sidechain.Anchor(query_mol, skeleton_contact, query_mcs_indices)
        all_sidechains.append(new_sidechain)

      
    

#- for each fragment define contact and anchor in terms of mol
#- order sidechains - reconstruction has to be done in order of distance from contact - get an order

#- steal xyz for mcs
#- measure internals on query side_chains (using contact+anchor based on mol indices)
#- reconstruct based on internals obtained in previous step - only the contact atom takes dihedrals from template

Chem.rdmolops.GetShortestPath(mol,id1,id2)
FindAtomEnvironmentOfRadiusN


mol1.HasSubstructMatch(my_mcs)

