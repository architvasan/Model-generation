"""This module contains functions docking a molecule to a receptor using Openeye.

The code is adapted from this repository: https://github.com/inspiremd/Model-generation
"""
import MDAnalysis as mda
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from functools import cache, partial
from pathlib import Path
from typing import List, Optional
import numpy as np
from openeye import oechem, oedocking, oeomega
import pandas as pd
from mpi4py import MPI

#from pipt.docking_utils import smi_to_structure
from utils import exception_handler
import os

def from_mol(mol, isomer=True, num_enantiomers=1):
    """
    Generates a set of conformers as an OEMol object
    Inputs:
        mol is an OEMol
        isomers is a boolean controling whether or not the various diasteriomers of a molecule are created
        num_enantiomers is the allowable number of enantiomers. For all, set to -1
    """
    # Turn off the GPU for omega
    omegaOpts = oeomega.OEOmegaOptions()
    omegaOpts.GetTorDriveOptions().SetUseGPU(False)
    omega = oeomega.OEOmega(omegaOpts)

    out_conf = []
    if not isomer:
        ret_code = omega.Build(mol)
        if ret_code == oeomega.OEOmegaReturnCode_Success:
            out_conf.append(mol)
        else:
            oechem.OEThrow.Warning(
                "%s: %s" % (mol.GetTitle(), oeomega.OEGetOmegaError(ret_code))
            )

    elif isomer:
        for enantiomer in oeomega.OEFlipper(mol.GetActive(), 12, True):
            enantiomer = oechem.OEMol(enantiomer)
            ret_code = omega.Build(enantiomer)
            if ret_code == oeomega.OEOmegaReturnCode_Success:
                out_conf.append(enantiomer)
                num_enantiomers -= 1
                if num_enantiomers == 0:
                    break
            else:
                oechem.OEThrow.Warning(
                    "%s: %s" % (mol.GetTitle(), oeomega.OEGetOmegaError(ret_code))
                )
    return out_conf


def from_string(smiles, isomer=True, num_enantiomers=1):
    """
    Generates an set of conformers from a SMILES string
    """
    mol = oechem.OEMol()
    if not oechem.OESmilesToMol(mol, smiles):
        raise ValueError(f"SMILES invalid for string {smiles}")
    else:
        return from_mol(mol, isomer, num_enantiomers)


def from_structure(structure_file: Path) -> oechem.OEMol:
    """
    Generates an set of conformers from a SMILES string
    """
    mol = oechem.OEMol()
    ifs = oechem.oemolistream()
    if not ifs.open(str(structure_file)):
        raise ValueError(f"Could not open structure file: {structure_file}")

    if structure_file.suffix == ".pdb":
        oechem.OEReadPDBFile(ifs, mol)
    elif structure_file.suffix == ".sdf":
        oechem.OEReadMDLFile(ifs, mol)
    else:
        raise ValueError(f"Invalid structure file extension: {structure_file}")

    return mol


def select_enantiomer(mol_list):
    return mol_list[0]


def dock_conf(receptor, mol, max_poses: int = 1):
    dock = oedocking.OEDock()
    dock.Initialize(receptor)
    lig = oechem.OEMol()
    #err = dock.DockMultiConformerMolecule(mol, max_poses)
    err = dock.DockMultiConformerMolecule(lig, mol, max_poses)
    #print(err)
    # the above line outputs the error
    return dock, lig

# Returns an array of length max_poses from above. This is the range of scores
def ligand_scores(dock, lig):
    return [dock.ScoreLigand(conf) for conf in lig.GetConfs()]


def best_dock_score(dock, lig):
    return ligand_scores(dock, lig)#[0]


def write_ligand(ligand, output_dir: Path, smiles: str) -> None:
    # TODO: If MAX_POSES != 1, we should select the top pose to save
    ofs = oechem.oemolostream()
    for it, conf in enumerate(list(ligand.GetConfs())):
    #conf = list(ligand.GetConfs())[0]
        if ofs.open(f'{str(output_dir)}/{smiles}/{it}.pdb'):
            oechem.OEWriteMolecule(ofs, conf)
            ofs.close()
    return
    raise ValueError(f"Could not write ligand to {output_path}")


def write_receptor(receptor, output_path: Path) -> None:
    ofs = oechem.oemolostream()
    #print(ofs)
    if ofs.open(str(output_path)):
        #print("Okie!")
        mol = oechem.OEMol()
        contents = receptor.GetComponents(mol)#Within
        #print("Yo!")
        #print(mol)
        oechem.OEWriteMolecule(ofs, mol)
        ofs.close()
    return
    raise ValueError(f"Could not write receptor to {output_path}")


#@cache  # Only read the receptor once
def read_receptor(receptor_oedu_file: Path):
    """Read the .oedu file into a GraphMol object."""
    receptor = oechem.OEDesignUnit()
    oechem.OEReadDesignUnit(str(receptor_oedu_file), receptor)
    return receptor

#@cache
def create_proteinuniv(protein_pdb):
    protein_universe = mda.Universe(protein_pdb)
    return protein_universe

def create_complex(protein_universe, ligand_pdb):
    u1 = protein_universe
    u2 = mda.Universe(ligand_pdb)
    u = mda.core.universe.Merge(u1.select_atoms("chainID A and not resname AKG NME HOH"), u2.atoms)#, u3.atoms)
    return u

def create_trajectory(protein_universe, ligand_dir, output_file_name):
    import MDAnalysis as mda
    ligand_files = sorted(os.listdir(ligand_dir))
    comb_univ_1 = create_complex(protein_universe, f'{ligand_dir}/{ligand_files[0]}').select_atoms("all")
    with mda.Writer(output_file_name, comb_univ_1.n_atoms,) as w:
        for it, ligand_file in enumerate(ligand_files):
            comb_univ = create_complex(protein_universe, f'{ligand_dir}/{ligand_file}') 
            w.write(comb_univ)    # write a whole universe
            os.remove(f'{ligand_dir}/{ligand_file}')
    return

#@exception_handler(default_return=0.0)
def run_docking(
    smiles: str, receptor_oedu_file: Path, max_confs: int, score_cutoff: float, temp_dir: Path, out_lig_dir: Path, temp_storage: Optional[Path] = None
) -> float:
    """Run OpenEye docking on a single ligand.

    Parameters
    ----------
    smiles : ste
        A single SMILES string.
    receptor_oedu_file : Path
        Path to the receptor .oedu file.
    max_confs : int
        Number of ligand poses to generate
    temp_lig_dir : Path
        Temporary directory to write individual ligand poses
    out_lig_path : Path
        Location to output ligand_protein top n poses --> single pdb
    out_rec_path : Path
        Where to write output receptor file
    temp_storage : Path
        Path to the temporary storage directory to write structures to,
        if None, use the current working Python's built in temp storage.

    Returns
    -------
    float
        The docking score of the best conformer.
    """

    try:
        print(smiles)
        conformers = select_enantiomer(from_string(smiles))

        # Read the receptor to dock to
        receptor = read_receptor(receptor_oedu_file)

        # Dock the ligand conformers to the receptor
        #print(receptor)
        #print(conformers)
        #print(max_confs)
        dock, lig = dock_conf(receptor, conformers, max_poses=max_confs)

        #print(receptor_oedu_file)
        # Get the docking scores
        best_score = best_dock_score(dock, lig)
        # If receptor pdb file doesnt exist then create one
        #out_rec_path = f'{temp_dir}/{os.path.splitext(os.path.basename(str(receptor_oedu_file))[0])}.pdb'
        out_rec_path = f'{temp_dir}/{Path(receptor_oedu_file).stem}.pdb'
        #print(out_rec_path)

        if not os.path.isfile(out_rec_path):
            write_receptor(receptor, out_rec_path)
        if best_score[0]<score_cutoff:
            if not os.path.isdir(f'{temp_dir}/{smiles}'):
                os.mkdir(f'{temp_dir}/{smiles}')
            write_ligand(lig, temp_dir, smiles)
            protein_universe = create_proteinuniv(f'{out_rec_path}')
            create_trajectory(protein_universe, f'{temp_dir}/{smiles}', f'{out_lig_dir}/{Path(receptor_oedu_file).stem}.{smiles}.pdb')
        return np.mean(best_score[0:10])
    except:
        return 0

def run_mpi_docking(
    smiles_list: List[str],
    receptor_oedu_file: Path,
    max_confs: int,
    score_cutoff: float,
    temp_dir: Path,
    out_lig_dir: Path,
    out_csv_pattern: str
    ) -> float:

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank() 
    size = comm.Get_size()
    
    smiles_list_new = np.array_split(smiles_list, size)[rank]
    scores = [run_docking(smi,receptor_oedu_file, max_confs, score_cutoff, temp_dir, out_lig_dir) for smi in smiles_list_new]
    df_new = pd.DataFrame({'SMILES':smiles_list_new,'Scores':scores})
    df_new.to_csv(f'{out_csv_pattern}.{rank}.csv', index=False)
    return 1

def _run_parallel_docking(
    smiles_batch: List[str],
    receptor_oedu_file: Path,
    max_confs: int,
    temp_dir: Path,
    out_lig_dir: Path,
    num_workers: int,
    temp_storage: Optional[Path] = None,
) -> List[float]:
    # Run the docking computation
    #run_docking(smiles, receptor_oedu_file, max_confs, temp_dir, out_lig_dir)
    #scores = run_docking(smiles, receptor_oedu_file, max_confs, out_lig_dir, out_lig_path)

    worker_fn = partial(
        run_docking, receptor_oedu_file=receptor_oedu_file, max_confs=max_confs, temp_dir=temp_dir, out_lig_dir=out_lig_dir, temp_storage=temp_storage
    )
    docking_scores = []
    with ProcessPoolExecutor(max_workers=num_workers) as pool:
        for score in pool.map(worker_fn, smiles_batch):
            docking_scores.append(score)
    return docking_scores


def _create_complex(protein_pdb, ligand_dir):
    import numpy as np
    import os
    
    # Load protein structure (single frame)
    protein = mda.Universe(protein_pdb)
    
    # Specify the mass of Mn
    mn_mass = 54.93804  # Replace with the actual mass of Mn
    
    # Assign mass to Mn atom
    mn_atom_indices = protein.select_atoms('type MN').indices
    protein.atoms[mn_atom_indices].masses = mn_mass

    # Create an output directory for writing the combined trajectories
    output_dir = 'output_combined_trajectories'
    os.makedirs(output_dir, exist_ok=True)
    
    # Loop through ligand PDB files
    ligand_files = sorted(os.listdir(ligand_dir))
    for ligand_file in ligand_files:
        if ligand_file.endswith('.pdb'):
            # Load ligand structure (single frame)
            #print(ligand_file)
            ligand = mda.Universe(os.path.join(ligand_dir, ligand_file))
    
            # Create an output PDB file for each combination of protein and ligand frame
            output_file = os.path.join(output_dir, f'combined_{ligand_file}')
            with mda.Writer(output_file, protein.atoms.n_atoms + ligand.atoms.n_atoms) as pdb_writer:
                # Concatenate protein and ligand coordinates
                combined_coordinates = np.concatenate([protein.atoms.positions, ligand.atoms.positions], axis=0)
                #mn_atom_indices = combined_coordinates.select_atoms('name MN').indices
                #combined_coordinates[mn_atom_indices].masses = mn_mass

                # Create a temporary Universe with the combined coordinates
                # Create a temporary AtomGroup with the combined coordinates
                temp_atomgroup = mda.AtomGroup(temp_universe.atoms, universe=protein_universe)

                # Update the coordinates of the temporary AtomGroup
                temp_atomgroup.positions = combined_coordinates

                # Write the frame to the output PDB file
                pdb_writer.write(temp_atomgroup)



                #frame = mda.Merge(protein, ligand)
                #frame.atoms.positions = combined_coordinates
                ## Write the frame to the output PDB file
                #pdb_writer.write(frame.atoms)

                #temp_universe = mda.Universe.empty(n_atoms=len(combined_coordinates))
                #temp_universe.atoms.positions = combined_coordinates
                #temp_universe.add_TopologyAttr('name', 'temp_group')
                #temp_universe.load_new(combined_coordinates, format='array')
                # Write the frame to the output PDB file
                #pdb_writer.write(temp_universe.atoms) 
                # Write the combined coordinates to the output PDB file
                #pdb_writer.write(combined_coordinates)





if False:
    if __name__ == "__main__":
        from argparse import ArgumentParser
    
        parser = ArgumentParser()
        parser.add_argument(
            "-r", "--receptor", type=Path, required=True, help="Receptor .oedu file"
        )
        parser.add_argument(
            "-s", "--smiles", type=Path, required=True, help="Ligand SMILES .smi file"
        )
        parser.add_argument(
            "-o", "--output", type=Path, required=True, help="Output .csv file"
        )
        parser.add_argument(
            "-c", "--confs", type=Path, required=True, help="Max conformations for ligand docking"
        )
        parser.add_argument(
            "-t", "--temp", type=Path, required=True, help="Temp dir storage ligand pdbs"
        )

        parser.add_argument(
            "-d", "--outdir", type=Path, required=True, help="Output prot-lig pose trajectories"
        )

        parser.add_argument(
            "-g", "--storage", type=Path, default=None
        )
        parser.add_argument(
            "-n", "--num_workers", type=int, default=1
        )
        args = parser.parse_args()

        #smiles_batch: List[str],
        #receptor_oedu_file: Path,
        #max_confs: int,
        #temp_dir: Path,
        #out_lig_dir: Path,
        #num_workers: int,
        #temp_storage: Optional[Path] = None,
   
        smiles_list = args.smiles.read_text().split("\n")
    
        docking_scores = run_parallel_docking(
            smiles_list, args.receptor, args.confs, args.temp, args.outdir, args.num_workers, args.storage
        )
    
        # Format the output file
        file_contents = "SMILES, DockingScore\n"
        file_contents += "\n".join(
            f"{smiles},{score}" for smiles, score in zip(smiles_list, docking_scores)
        )
        file_contents += "\n"
    
        # Write the output file
        with open(args.output, "w") as f:
            f.write(file_contents)


df_smiles = pd.read_csv('input/6ie2_MNbound_top100.smi')['SMILES']
recep_file = 'input/6ie2_lig_akg.oedu'
max_confs = 100
score_cutoff = 0
temp_dir = 'lig_confs'
output_traj = 'output_combined_trajectories'
score_pattern = '6ie2_scores'

run_mpi_docking(df_smiles, recep_file, max_confs, score_cutoff, temp_dir, output_traj, score_pattern)

