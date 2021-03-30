import sys
import glob

# pyrosetta
import pyrosetta
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta.utility import vector1_std_string as vec1_string
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector

try:
    from pyRIF import RotamerInteractionField
except:
    sys.path.append('../pyRIF')
    from pyRIF import RotamerInteractionField

# init pyrosetta with whatever you want
pyrosetta.init(
    '-mute all '
    '-beta '
    '-in:file:extra_res_fa /home/psalveso/git_repositories/peptim_stapler/peptim_stapler/_inputs/CYB.params '
    '/home/psalveso/git_repositories/peptim_stapler/peptim_stapler/_inputs/NMX.params '
    '/home/psalveso/git_repositories/peptim_stapler/peptim_stapler/_inputs/B3Z.params '
    '/home/psalveso/projects/BRD/params/CYY.params '
    '/home/psalveso/projects/BRD/params/KDD.params '
)

# location of RIF files
L_AA_RIF = {
    'HDF5'   : '/net/scratch/psalveso/BRDs/RIF/6u8m/rifgen/py_rif.h5',
    'rots'   : '/net/scratch/psalveso/BRDs/RIF/6u8m/rifgen/rotamer_index_spec.txt',
    'target' : '/net/scratch/psalveso/BRDs/RIF/6u8m/rifgen/rif_64_6u8m_copy1_sca0.8_noKR.rif.gz_target.pdb.gz'
}

D_AA_RIF = {
    'HDF5'   : '/net/scratch/psalveso/BRDs/RIF/6u8m/rifgen_d/py_rif.h5',
    'rots'   : '/net/scratch/psalveso/BRDs/RIF/6u8m/rifgen_d/rotamer_index_spec.txt',
    'target' : '/net/scratch/psalveso/BRDs/RIF/6u8m/rifgen_d/rif_64_6u8m_copy1_sca0.8_noKR.rif.gz_target.pdb.gz',
}

# some test PDBs
PDBs = glob.glob(f'/home/psalveso/projects/BRD/peptim_cyclize_input_complex/6u8m/PEP/*.pdb')

# setup residue selector to select the target chain
chain_string = vec1_string()
chain_string.append('1')
target_selector = ChainSelector()
target_selector.set_chain_strings(chain_string)

# setup residue selector to select the binder chain
chain_string = vec1_string()
chain_string.append('2')
binder_selector = ChainSelector()
binder_selector.set_chain_strings(chain_string)


# make the RIF
RIF = RotamerInteractionField(
    L_AA_RIF_kwargs=L_AA_RIF,
    D_AA_RIF_kwargs=D_AA_RIF,
    residue_selector=binder_selector,
    target_selector=target_selector,
    min_RIF_hits=1,
    min_RIF_score=-0.5,
)
print('L/D RIF')
for PDB in PDBs:
    pose = pose_from_pdb(PDB)
    # apply the RIF
    STATUS, RIF_score, sequence_mapping = RIF.apply(pose)
    # test if pose hit enough RIF interactions, with low enough score
    if STATUS:
        print(f'rif score: {RIF_score}\n{sequence_mapping}')
        # do the remainder of your protocl in here....
        #
        # if you want the RIF residues
        # mutate them onto pose using sequence mapping


# RIF just with canonical L-AAs
RIF = RotamerInteractionField(
    L_AA_RIF_kwargs=L_AA_RIF,
    residue_selector=binder_selector,
    target_selector=target_selector,
)
print('L RIF')
for PDB in PDBs:
    pose = pose_from_pdb(PDB)
    # apply the RIF
    STATUS, RIF_score, sequence_mapping = RIF.apply(pose)
    # test if pose hit enough RIF interactions, with low enough score
    if STATUS: print(f'rif score: {RIF_score}\n{sequence_mapping}')


# RIF just with canonical D-AAs
RIF = RotamerInteractionField(
    D_AA_RIF_kwargs=D_AA_RIF,
    residue_selector=binder_selector,
    target_selector=target_selector,
)
print('D RIF')
for PDB in PDBs:
    pose = pose_from_pdb(PDB)
    # apply the RIF
    STATUS, RIF_score, sequence_mapping = RIF.apply(pose)
    # test if pose hit enough RIF interactions, with low enough score
    if STATUS: print(f'rif score: {RIF_score}\n{sequence_mapping}')

# RIF with hotspot searching
# make the RIF
RIF = RotamerInteractionField(
    L_AA_RIF_kwargs=L_AA_RIF,
    D_AA_RIF_kwargs=D_AA_RIF,
    residue_selector=binder_selector,
    target_selector=target_selector,
    min_HOTSPOT_hits=1,
    HOTSPOTs_in_RIF=8,
)
print('L/D Hotspot RIF')
for PDB in PDBs:
    pose = pose_from_pdb(PDB)
    # apply the RIF
    STATUS, RIF_score, sequence_mapping = RIF.apply(pose)
    # test if pose hit enough RIF interactions, with low enough score
    if STATUS:
        print(f'rif score: {RIF_score}\n{sequence_mapping}')
        # do the remainder of your protocl in here....
        #
        # if you want the RIF residues
        # mutate them onto pose using sequence mapping
