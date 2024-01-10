# Workflow 0
from openeye import oechem
import pandas as pd
import numpy as np
from impress_md import interface_functions
import argparse
import os
from mpi4py import MPI
from tqdm import tqdm

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
world_size = comm.Get_size()

import signal

WORKTAG, DIETAG = 11, 13
CHUNKSIZE=5

class timeout:
    def __init__(self, seconds=1, error_message='Timeout'):
        self.seconds = seconds
        self.error_message = error_message
    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)
    def __exit__(self, type, value, traceback):
        signal.alarm(0)

def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help='input csv for smiles', required=True, type=str)
    parser.add_argument("-o", help='output file for data', required=True, type=str)
    parser.add_argument("-r", help='receptor file', required=True, type=str)
    parser.add_argument('-v', help='verbose', action='store_true')
    parser.add_argument("-n", type=int, default=1)
    parser.add_argument("-l", type=str, default=None, required=False)
    return parser.parse_args()


def get_root_protein_name(file_name):
    return file_name.split("/")[-1].split(".")[0]


def get_smiles_col(col_names):
    return int(np.where(['smile' in s.lower() for s in col_names])[0][0])


def get_ligand_name_col(col_names):
    return int(np.where(['id' in s.lower() or 'title' in s.lower() or "name" in s.lower() for s in col_names])[0][0])


def master():
    smiles_file = pd.read_csv(input_smiles_file)
    smiles_file.sample(frac=1)
    columns = smiles_file.columns.tolist()
    smiles_col = get_smiles_col(columns)
    name_col = get_ligand_name_col(columns)
    it = 0 
    for pos in range(0, smiles_file.shape[0], CHUNKSIZE):
        print(it)
        it+=1
        poss = list(range(pos, pos + CHUNKSIZE))
        smiles = smiles_file.iloc[pos: min(pos + CHUNKSIZE, smiles_file.shape[0]), smiles_col]
        ligand_name = smiles_file.iloc[pos: min(pos + CHUNKSIZE, smiles_file.shape[0]), name_col]
        status = MPI.Status()
        comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        rank_from = status.Get_source()
        data = zip(poss, smiles, ligand_name)
        comm.send(data, dest=rank_from, tag=23)

def slave():
    data_dock = {}
    data_dock['smiles'] = []
    data_dock['score'] = []
    docker, receptor = interface_functions.get_receptor(target_file, use_hybrid=use_hybrid,
                                                        high_resolution=high_resolution)

    comm.send([], dest=0, tag=11)
    poss = comm.recv(source=0, tag=MPI.ANY_TAG)

    while True:

        for (pos, smiles, ligand_name) in poss:
            try:
                with timeout(seconds=180):
                    score, res, ligand = interface_functions.RunDocking_(smiles,
                                                                         dock_obj=docker,
                                                                         pos=pos,
                                                                         name=ligand_name,
                                                                         target_name=pdb_name,
                                                                         force_flipper=force_flipper)
                    #data_dock['smiles'].append(smiles)
                    #data_dock['score'].append(score)

                    with open(f"all_docking/docking_{rank}.dat", "a") as myfile:
                        myfile.write(f"{smiles}, {score}\n")

                    if args.v:
                        print("RANK {}:".format(rank), res, end='')

                    #if ofs and ligand is not None:
                    #    oechem.OEWriteMolecule(ofs, ligand)
            except TimeoutError:
                print("TIMEOUT", smiles, ligand_name)
                continue
        comm.send([], dest=0, tag=11)

        status = MPI.Status()

        poss = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        if status.Get_tag() == DIETAG:
            break

    #if ofs is not None:
    #    ofs.close()
    comm.Barrier()


if __name__ == '__main__':
    args = getargs()
    input_smiles_file = args.i
    target_file = args.r  # twenty of these
    basename, file_ending = ".".join(args.o.split(".")[:-1]), args.o.split(".")[-1]
    output_poses = basename + "_" + str(rank) + "." + file_ending
    pdb_name = get_root_protein_name(target_file)
    ## setting don't change
    use_hybrid = True
    force_flipper = False
    high_resolution = True
    #ofs = oechem.oemolostream(output_poses)

    with open(f"all_docking/docking_{rank}.dat", "w") as myfile:
        myfile.write("smiles,score\n")

    if rank == 0:
        master()
    else:
        slave()
