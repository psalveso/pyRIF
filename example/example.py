import sys
import glob
sys.path.append('/home/psalveso/git_repositories/pyRIF/pyRIF')

# pyrosetta
import pyrosetta
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta.utility import vector1_std_string as vec1_string
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector


from pyRIF import RotamerInteractionField

# init pyrosetta with whatever you want
pyrosetta.init(
    '-mute all '
    '-in:file:extra_res_fa /home/psalveso/git_repositories/peptim_stapler/peptim_stapler/_inputs/CYB.params '
    '/home/psalveso/git_repositories/peptim_stapler/peptim_stapler/_inputs/NMX.params '
    '/home/psalveso/git_repositories/peptim_stapler/peptim_stapler/_inputs/B3Z.params '
    #'/home/psalveso/git_repositories/peptim_stapler/peptim_stapler/_inputs/401.params '
)

# location of RIF files
L_AA_RIF = {
    'HDF5'   : '/net/scratch/psalveso/peptim_cyclize_stapler_outputs/TEAD/RIF_HOTSPOTS/site_65-69/rifdock/py_rif.h5',
    'rots'   : '/net/scratch/psalveso/peptim_cyclize_stapler_outputs/TEAD/RIF_HOTSPOTS/site_65-69/rifgen/rotamer_index_spec.txt',
    'target' : '/net/scratch/psalveso/peptim_cyclize_stapler_outputs/TEAD/RIF_HOTSPOTS/site_65-69/rifgen/rif_64_TEAD_sca0.8_noKR.rif.gz_target.pdb.gz'
}

D_AA_RIF = {
    'HDF5'   : '/net/scratch/psalveso/peptim_cyclize_stapler_outputs/RIF/rifdock_TEAD_D/py_rif.h5',
    'rots'   : '/net/scratch/psalveso/peptim_cyclize_stapler_outputs/RIF/rifgen_TEAD_D/rotamer_index_spec.txt',
    'target' : '/net/scratch/psalveso/peptim_cyclize_stapler_outputs/RIF/rifgen_TEAD_D/rif_64_TEAD_sca0.8_noKR.rif.gz_target.pdb.gz',
}


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
)


PDBs = glob.glob(f'/net/scratch/psalveso/peptim_cyclize_stapler_outputs/TEAD/test_RIF_sacffolds/*.pdb')

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
for PDB in PDBs:
    pose = pose_from_pdb(PDB)
    # apply the RIF
    STATUS, RIF_score, sequence_mapping = RIF.apply(pose)
    # test if pose hit enough RIF interactions, with low enough score
    if STATUS: print(f'rif score: {RIF_score}\n{sequence_mapping}')
