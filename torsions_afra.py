from __future__ import print_function
import mdtraj as md
import numpy as np
import sys
import subprocess
import math
import mdtraj as md
import numpy as np
import math
import pickle
from Bio.PDB import PDBParser
from Bio import PDB
import pandas as pd
import itertools
from collections import OrderedDict

nres = []
resname = []
ang = []
i_list = []


def file0(traj):

    with open("list_structures.temp") as file:
        lines = file.readlines()
    for i in lines:
        parser = PDBParser(PERMISSIVE=1)
        traj = i.split()[1]
        structure = parser.get_structure("lig", traj)

        residues = [r for r in structure.get_residues()]

        for i in range(len(residues)):
            i = i+1
            nres.append(i)
        for i in residues:
            name = i.get_resname()
            resname.append(name)
            # print(name)

        for each in itertools.combinations(residues, 2):
            seq1 = each[0].get_id()[1]
            seq2 = each[1].get_id()[1]
            seq1_name = each[0].get_resname()
            seq2_name = each[1].get_resname()

            if (seq1_name == 'U' or seq1_name == 'C'):
                Co1 = each[0]["N1"].get_vector()
                Co2 = each[0]["O4'"].get_vector()
                Co3 = each[0]["C3'"].get_vector()
                Co4 = each[0]["O2'"].get_vector()
                coord = [Co1, Co2, Co3, Co4]
                print(coord)

            else:
                Co1 = each[0]["N9"].get_vector()
                Co2 = each[0]["O4'"].get_vector()
                Co3 = each[0]["C3'"].get_vector()
                Co4 = each[0]["O2'"].get_vector()
                coord = [Co1, Co2, Co3, Co4]
                print(coord)
            if (seq2-1) == seq1:
                angle = PDB.calc_dihedral(Co1, Co2, Co3, Co4)
            ang.append(angle)
            print(angle)
        ang2 = list(OrderedDict.fromkeys(ang))
        ang3 = pd.DataFrame(ang2)
        angle_pi = ang3.apply(lambda x: (x*math.pi)*180 /
                              math.pi if x[0] < 0 else 0, axis=1)
    np.savetxt('dih_NbO4pC3pO2p_trj_test.dat', angle_pi, fmt='%1.2f')
