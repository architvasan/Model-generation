import sys, os
import numpy as np
from contextlib import contextmanager

@contextmanager
def working_directory(directory):
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)

def RunDocking(smiles, inpath, outpath, padding=4):
    from . import conf_gen
    from . import dock_conf
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    confs = conf_gen.SelectEnantiomer(conf_gen.FromString(smiles))
    # This receptor can be pre-compiled to an oeb. It speeds things up
    filename, file_extension = os.path.splitext(inpath)
    #if file_extension == ".oeb":
    #    receptor = dock_conf.PrepareReceptorFromBinary(inpath)
    #else: # else it is a pdb
    #    receptor = dock_conf.PrepareReceptor(inpath,padding,outpath)
    
    dock, lig, receptor = dock_conf.DockConf(inpath,file_extension,confs[0],MAX_POSES=1)
    
    # Currently we generate 200 conformers for each ligand, but only take
    #   the best pose, as scored by Openeye. It may be useful to consider
    #   something about the range of poses.

    dock_conf.WriteStructures(receptor, lig, f'{outpath}/apo.pdb', f'{outpath}/lig.pdb')
    with open(f'{outpath}/metrics.csv','w+') as metrics:
        metrics.write("Dock,Dock_U\n")
        metrics.write("{},{}\n".format(dock_conf.BestDockScore(dock,lig),0))
    # # If you uncomment the three lines below, it will save an image of the 2D
    #   molecule. This is useful as a sanity check.
    # from openeye import oedepict
    # oedepict.OEPrepareDepiction(lig)
    # oedepict.OERenderMolecule(f'{outpath}/lig.png',lig)


def RunDocking_(smiles, inpath, outpath, padding=4):
    from . import conf_gen
    from . import dock_conf
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    confs = conf_gen.SelectEnantiomer(conf_gen.FromString(smiles))
    # This receptor can be pre-compiled to an oeb. It speeds things up
    filename, file_extension = os.path.splitext(inpath)
    #if file_extension == ".oeb":
    #    receptor = dock_conf.PrepareReceptorFromBinary(inpath)
    #else: # else it is a pdb
    #    receptor = dock_conf.PrepareReceptor(inpath,padding,outpath)

    dock, lig, receptor = dock_conf.DockConf("input/receptor.oeb",confs,MAX_POSES=1)

    # Currently we generate 200 conformers for each ligand, but only take
    #   the best pose, as scored by Openeye. It may be useful to consider
    #   something about the range of poses.

    with open(f'{outpath}/metrics.csv','w+') as metrics:
        metrics.write("Dock,Dock_U\n")
        metrics.write("{},{}\n".format(dock_conf.BestDockScore(dock,lig),0))
    dock_conf.WriteStructures(receptor, lig, f'{outpath}/apo.pdb', f'{outpath}/lig.pdb')
    # # If you uncomment the three lines below, it will save an image of the 2D
    #   molecule. This is useful as a sanity check.
    # from openeye import oedepict
    # oedepict.OEPrepareDepiction(lig)
    # oedepict.OERenderMolecule(f'{outpath}/lig.png',lig)
    return dock_conf.BestDockScore(dock,lig)

def ParameterizeOE(path):
    """
    Reads in the PDB from 'RunDocking' and outputs 'charged.mol2' of the ligand
    Then runs antechamber to convert this to coordinate (.inpcrd) and 
    parameter (.prmtop) files.
    """
    from openeye import oechem, oeomega, oequacpac
    mol = oechem.OEMol()
    ifs = oechem.oemolistream()
    if ifs.open(f'{path}/lig.pdb'):
        oechem.OEReadMolecule(ifs,mol)
        ifs.close()
    if not oequacpac.OEAssignCharges(mol,oequacpac.OEAM1BCCCharges()):
        raise(RuntimeError("OEAssignCharges failed."))
    ofs = oechem.oemolostream()
    if ofs.open(f'{path}/charged.mol2'):
        oechem.OEWriteMolecule(ofs,mol)
    
    import subprocess
    with working_directory(path):
        subprocess.check_output(f'antechamber -i lig.pdb -fi pdb -o lig.mol2 -fo mol2 -pf y -an y -a charged.mol2 -fa mol2 -ao crg',shell=True)
        subprocess.check_output(f'parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod',shell=True)
        # Wrap tleap
        with open(f'leap.in','w+') as leap:
            leap.write("source leaprc.protein.ff14SBonlysc\n")
            leap.write("source leaprc.gaff\n")
            leap.write("source leaprc.water.tip3p\n")
            leap.write("set default PBRadii mbondi3\n")
            leap.write("rec = loadPDB apo.pdb # May need full filepath?\n")
            leap.write("saveAmberParm rec apo.prmtop apo.inpcrd\n")
            leap.write("lig = loadmol2 lig.mol2\n")
            leap.write("loadAmberParams lig.frcmod\n")
            leap.write("com = combine {rec lig}\n")
            leap.write("saveAmberParm lig lig.prmtop lig.inpcrd\n")
            leap.write("saveAmberParm com com.prmtop com.inpcrd\n")
            leap.write("savepdb com com.pdb\n")
            leap.write("addIons2 rec NA 0\n")
            leap.write("addIons2 rec CL 0\n")
            leap.write("solvateBox rec TIP3PBOX 14\n")
            leap.write("savepdb rec apo_sol.pdb\n")
            leap.write("saveAmberParm rec apo_sol.prmtop apo_sol.inpcrd\n")
            leap.write("addIons2 lig NA 0\n")
            leap.write("addIons2 lig CL 0\n")
            leap.write("solvateBox lig TIP3PBOX 14\n")
            leap.write("savepdb lig lig_sol.pdb\n")
            leap.write("saveAmberParm lig lig_sol.prmtop lig_sol.inpcrd\n")
            leap.write("addIons2 com NA 0\n")
            leap.write("addIons2 com CL 0\n")
            leap.write("solvateBox com TIP3PBOX 14\n")
            leap.write("savepdb com com_sol.pdb\n")
            leap.write("saveAmberParm com com_sol.prmtop com_sol.inpcrd\n")
            leap.write("quit\n")
        subprocess.check_output(f'tleap -f leap.in',shell=True)


def ParameterizeAMBER(path):
    """
    Alternative method for parameterizing the system. It uses sqm and runs much slower than
    ParameterizeOE, which works through OpenEye.
    This function is pretty much a wrapper for antechamber & tleap. 
    I've kept it as a backup in case there are issues with OpenEye. 
    """
    import subprocess
    with working_directory(path):
        subprocess.check_output(f'antechamber -i lig.pdb -fi pdb -o lig.mol2 -fo mol2 -c bcc -pf y -an y',shell=True)
        subprocess.check_output(f'parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod',shell=True)
        with open(f'leap.in','w+') as leap:
            leap.write("source leaprc.protein.ff14SBonlysc\n")
            leap.write("source leaprc.gaff\n")
            leap.write("source leaprc.water.tip3p\n")
            leap.write("set default PBRadii mbondi3\n")
            leap.write("rec = loadPDB apo.pdb # May need full filepath?\n")
            leap.write("saveAmberParm rec apo.prmtop apo.inpcrd\n")
            leap.write("lig = loadmol2 lig.mol2\n")
            leap.write("loadAmberParams lig.frcmod\n")
            leap.write("com = combine {rec lig}\n")
            leap.write("saveAmberParm lig lig.prmtop lig.inpcrd\n")
            leap.write("saveAmberParm com com.prmtop com.inpcrd\n")
            leap.write("savepdb com com.pdb\n")
            leap.write("addIons2 rec NA 0\n")
            leap.write("addIons2 rec CL 0\n")
            leap.write("solvateBox rec TIP3PBOX 14\n")
            leap.write("savepdb rec apo_sol.pdb\n")
            leap.write("saveAmberParm rec apo_sol.prmtop apo_sol.inpcrd\n")
            leap.write("addIons2 lig NA 0\n")
            leap.write("addIons2 lig CL 0\n")
            leap.write("solvateBox lig TIP3PBOX 14\n")
            leap.write("savepdb lig lig_sol.pdb\n")
            leap.write("saveAmberParm lig lig_sol.prmtop lig_sol.inpcrd\n")
            leap.write("addIons2 com NA 0\n")
            leap.write("addIons2 com CL 0\n")
            leap.write("solvateBox com TIP3PBOX 14\n")
            leap.write("savepdb com com_sol.pdb\n")
            leap.write("saveAmberParm com com_sol.prmtop com_sol.inpcrd\n")
            leap.write("quit\n")
        subprocess.check_output(f'tleap -f leap.in',shell=True)
    

def RunMinimization(build_path, outpath, one_traj=False):
    """
    We are minimizing all three structures, then checking the potential energy using GB forcefields
    We could, alternatively, minimize the docked structure and then extract trajectories (1 frame long),
    more like a 1-trajectory mmgbsa.
    output is path used in "RunDocking". It has a metric.csv file.
    """
    from . import minimize
    success = True
    try:
        rec_energy = minimize.MinimizedEnergy(f'{build_path}/apo')
        lig_energy = minimize.MinimizedEnergy(f'{build_path}/lig')
        com_energy = minimize.MinimizedEnergy(f'{build_path}/com')
        diff_energy = com_energy - lig_energy - rec_energy
    except:
        success = False
    
    if one_traj:
        print("1-traj calculation not ready")
    # TODO: We could decide to do 1-trajectory mmgbsa. It would run about twice as fast as the
    #       current method. I think it would be less accurate, but maybe not. Look into the 1-traj
    #       method from Coveney papers if you want to implement this.

    with open(f'{outpath}/metrics.csv','r') as metrics:
        dat = metrics.readlines()
    with open(f'{outpath}/metrics.csv','w') as metrics:
        metrics.write(dat[0].replace('\n',',Minimize,Minimize_U\n'))
        if success:
            metrics.write(dat[1].replace('\n',',{},{}\n'.format(diff_energy,0)))
        else:
            metrics.write(dat[1].replace('\n',',NA,NA\n'))
def RunMinimization_(build_path, outpath, one_traj=False):
    from . import minimize
    success = True
    try:
        rec_energy = minimize.MinimizedEnergy(f'{build_path}/apo')
        lig_energy = minimize.MinimizedEnergy(f'{build_path}/lig')
        com_energy = minimize.MinimizedEnergy(f'{build_path}/com')
        diff_energy = com_energy - lig_energy - rec_energy
    except:
        success = False

    if one_traj:
        print("1-traj calculation not ready")
    # TODO: We could decide to do 1-trajectory mmgbsa. It would run about twice as fast as the
    #       current method. I think it would be less accurate, but maybe not. Look into the 1-traj
    #       method from Coveney papers if you want to implement this.
    with open(f'{outpath}/metrics.csv','r') as metrics:
        dat = metrics.readlines()
    with open(f'{outpath}/metrics.csv','w') as metrics:
        metrics.write(dat[0].replace('\n',',Minimize,Minimize_U\n'))
        if success:
            metrics.write(dat[1].replace('\n',',{},{}\n'.format(diff_energy,0)))
        else:
            metrics.write(dat[1].replace('\n',',NA,NA\n'))
    if success:
        return diff_energy
    else:
        return np.nan


def Simulation_explicit(inpath, outpath, nsteps, comp):
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    from . import minimize
    success = True
    try:
        potential = minimize.simulation(f'{inpath}/{comp}', f'{outpath}/{comp}', nsteps)
    except:
        success = False

    with open(f'{inpath}/metrics.csv','r') as metrics:
        dat = metrics.readlines()
    with open(f'{inpath}/metrics.csv','w') as metrics:
        metrics.write(dat[0].replace('\n',f',{comp}_U_mean_explicit\n'))
        if success:
            metrics.write(dat[1].replace('\n',',{}\n'.format(potential)))
        else:
            metrics.write(dat[1].replace('\n',',NA\n'))
    if success:
        return potential
    else:
        return np.nan

def Simulation_ESMACS(inpath, outpath, nsteps, replica, comp):
    inpath = os.path.join(inpath,comp)
    outpath = os.path.join(outpath,comp,'rep'+str(replica)) 
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    from . import minimize
    try:
        minimize.simulation_ESMACS(inpath, outpath, nsteps)
    except:
        print('\nrep{} simulation did not finish successfully!\n'.format(replica))

def RunMMGBSA_explicit(inpath, outpath, three_traj):
    """
    Extracts conformations from the trajectories and calculates potential energy using GB forcefields (MMGBSA energy).
    """
    from . import minimize
    import traj_parser as tp
    import os

    tmppath = outpath + '/tmp'
    if not os.path.exists(tmppath):
        os.makedirs(tmppath)

    if three_traj == 'yes':
        comp = 'com'
        tp.extractconfs(outpath, tmppath, comp)
        comp = 'apo'
        tp.extractconfs(outpath, tmppath, comp)
        comp = 'lig'
        tp.extractconfs(outpath, tmppath, comp)
    else:
        comp = 'com'
        tp.extractconfs(outpath, tmppath, comp)
        import glob
        for file in glob.glob(os.path.join(tmppath, 'com_MMGBSA_snap*.pdb')):
            name = file.split('/')[-1]
            rec_name = tmppath + '/' + name.replace('com', 'apo', 1)
            lig_name = tmppath + '/' + name.replace('com', 'lig', 1)
            ## This may not be the best way to split pdb but it works assuming that we only have one
            ## ligand which is also the last residue in the complex pdb. MDAnalysis can do this easily and in a
            ## robust manner, but I wanted to avoid dependencies, so coded it this way.
            tp.splitpdb(file, rec_name, lig_name)

    success = True
    try:
        comp = 'com'
        com_energy = minimize.MMGBSA(inpath, tmppath, comp)
        comp = 'apo'
        rec_energy = minimize.MMGBSA(inpath, tmppath, comp)
        comp = 'lig'
        lig_energy = minimize.MMGBSA(inpath, tmppath, comp)
        diff_energy = com_energy - lig_energy - rec_energy
    except:
        success = False

    with open(f'{inpath}/metrics.csv','r') as metrics:
        dat = metrics.readlines()
    with open(f'{inpath}/metrics.csv','w') as metrics:
        metrics.write(dat[0].replace('\n',',del_G_MMGBSA\n'))
        if success:
            metrics.write(dat[1].replace('\n',',{}\n'.format(diff_energy)))
        else:
            metrics.write(dat[1].replace('\n',',NA\n'))

    if success:
        return diff_energy
    else:
        return np.nan


def RunMMGBSA(inpath, outpath, niter=1000):
    """
    1 'iteration' corresponds to 1 ps.
    """
    from . import mmgbsa
    crds = {'lig':f'{inpath}/lig.inpcrd','apo':f'{inpath}/apo.inpcrd','com':f'{inpath}/com.inpcrd'}
    prms = {'lig':f'{inpath}/lig.prmtop','apo':f'{inpath}/apo.prmtop','com':f'{inpath}/com.prmtop'}
    
    enthalpies = mmgbsa.simulate(crds, prms, niter)
    # enthalpies is a list of energies from each iteration
    mmgbsa.subsample(enthalpies)
    # We subsample the enthalpies using a method from John Chodera that determines the equilibration
    #   and autocorrelation times. This allows us to extract an uncertainty.
    #   See the file mmgbsa.py or his package 'pymbar' for more detail.
    energies = mmgbsa.mmgbsa(enthalpies)
    
    with open(f'{outpath}/metrics.csv','r') as metrics:
        dat = metrics.readlines()
    with open(f'{outpath}/metrics.csv','w') as metrics:
        metrics.write(dat[0].replace('\n',',mmgbsa,mmgbsa_U\n'))
        metrics.write(dat[1].replace('\n',',{},{}\n'.format(energies[0]['diff'],energies[1]['diff'])))
    return energies


def RunMMGBSA_(inpath, outpath, niter=1000):
    """
    1 'iteration' corresponds to 1 ps.
    """
    from . import mmgbsa
    crds = {'lig':f'{inpath}/lig.inpcrd','apo':f'{inpath}/apo.inpcrd','com':f'{inpath}/com.inpcrd'}
    prms = {'lig':f'{inpath}/lig.prmtop','apo':f'{inpath}/apo.prmtop','com':f'{inpath}/com.prmtop'}

    enthalpies = mmgbsa.simulate(crds, prms, niter)
    # enthalpies is a list of energies from each iteration
    mmgbsa.subsample(enthalpies)
    # We subsample the enthalpies using a method from John Chodera that determines the equilibration
    #   and autocorrelation times. This allows us to extract an uncertainty.
    #   See the file mmgbsa.py or his package 'pymbar' for more detail.
    energies = mmgbsa.mmgbsa(enthalpies)

    with open(f'{outpath}/metrics.csv','r') as metrics:
        dat = metrics.readlines()
    with open(f'{outpath}/metrics.csv','w') as metrics:
        metrics.write(dat[0].replace('\n',',mmgbsa,mmgbsa_U\n'))
        metrics.write(dat[1].replace('\n',',{},{}\n'.format(energies[0]['diff'],energies[1]['diff'])))
    return energies[0]['diff']

def RunAlchemy(path, niter=2500, nsteps_per_iter=1000, nlambda=11):
    """
    Default is a 5 ns simulation with sampling every 2 ps
    """
    from . import alchemy
    [energy, err] = alchemy.SimulateAlchemy(path, niter, nsteps_per_iter, nlambda)
    with open(f'{path}/metrics.csv','r') as metrics:
        dat = metrics.readlines()
    with open(f'{path}/metrics.csv','w') as metrics:
        metrics.write(dat[0].replace('\n',',alchemy,alchemy_U\n'))
        metrics.write(dat[1].replace('\n',',{},{}\n'.format(energy,err)))
    return energy, err
