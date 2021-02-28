#! /usr/bin/env python3
# this is a python script for generating crdchg file from a pdb file saved from VMD.
# or generating crdchg file from a xyz file.
# Author: Wen-Kai Chen
# Email: wenkaichen@mail.bnu.edu.cn

# let's start!



def from_pdb_to_crdchg(filename, nmol, natoms, Azo_chg, wat_chg):
    with open(filename, 'r') as f:
        lines = f.readlines()
    resnum_old = 0
    snap_shot_index = 1
    snap_shot = "{nmol}  {natoms}\n-----\n".format(nmol=nmol,natoms=natoms)
    for line in lines[1:]:
        if line != 'END\n': # not reach to the end of one cluster.
            resnum_new = int(line.split()[5])
            if resnum_new != resnum_old: # this is a new residue

                if resnum_old > resnum_new or resnum_old == 0: # very first residue LIG/acetone
                    snap_shot += "0 1 1\n"
                    LIG = True
                    data=line.split()
                    mol = "{ele:4s}{x:15.8f}{y:15.8f}{z:15.8f}{chg:11.5f}\n".format(ele=data[2][0],
                                                                                 x=float(data[6]),
                                                                                 y=float(data[7]),
                                                                                 z=float(data[8]),
                                                                                 chg=Azo_chg[int(data[1])-1])

                else: # not very first residue WAT
                    snap_shot += "0 1 2\n"
                    LIG = False
                    data = line.split()
                    mol = "{ele:4s}{x:15.8f}{y:15.8f}{z:15.8f}{chg:11.5f}\n".format(ele=data[2][0],
                                                                                 x=float(data[6]),
                                                                                 y=float(data[7]),
                                                                                 z=float(data[8]),
                                                                                 chg=wat_chg[0])
                resnum_old = resnum_new

            else:
                if LIG:
                    data=line.split()
                    mol = "{ele:4s}{x:15.8f}{y:15.8f}{z:15.8f}{chg:11.5f}\n".format(ele=data[2][0],
                                                                                 x=float(data[6]),
                                                                                 y=float(data[7]),
                                                                                 z=float(data[8]),
                                                                                 chg=Azo_chg[int(data[1])-1])
                else:
                    data = line.split()
                    mol = "{ele:4s}{x:15.8f}{y:15.8f}{z:15.8f}{chg:11.5f}\n".format(ele=data[2][0],
                                                                                 x=float(data[6]),
                                                                                 y=float(data[7]),
                                                                                 z=float(data[8]),
                                                                                 chg=wat_chg[1])
            snap_shot += mol

        else:
            if resnum_old != 0: # not the very first snapshot
                with open("cluster_{}.crdchg".format(snap_shot_index), 'w') as f:
                    f.write(snap_shot)
                    snap_shot = "{nmol}  {natoms}\n-----\n".format(nmol=nmol,natoms=natoms)
                resnum_old = resnum_new
                snap_shot_index += 1
            else: # very first snapshot
                pass


def from_txt_to_crdchg( nmol, natoms, mol_list,
                       dic_mol_natom,dic_mol_element,dic_mol_chg,filename='coord.txt'):
    with open(filename, 'r') as f:
        lines = f.readlines()
    natom_ace = len(Azo_chg)
    natom_wat = len(wat_chg)
    current_mol = 0
    mol_num = []
    snap_shot = "{nmol}  {natoms}\n-----\n".format(nmol=nmol,natoms=natoms)
    mol_num = [0]
    for _,i in enumerate(mol_list):
        if _ == 0:
            mol_num.append(dic_mol_natom[i])
        else:
            mol_num.append(dic_mol_natom[i]+mol_num[_])
    for index, line in enumerate(lines):
        while True:
            if index < mol_num[current_mol+1]: # this is the correct molecular.
                atom_index = index - mol_num[current_mol]
                break
            else:
                current_mol += 1
        mol_name = mol_list[current_mol]
        if atom_index == 0:
            snap_shot += '0 1 {index}\n'.format(index=dic_mol_index[mol_name])
        snap_shot += "{ele:4s}{x:15.8f}{y:15.8f}{z:15.8f}{chg:11.5f}\n".format(ele=dic_mol_element[mol_name][atom_index],
                                                                               x=float(line.split()[0]),
                                                                               y=float(line.split()[1]),
                                                                               z=float(line.split()[2]),
                                                                               chg=dic_mol_chg[mol_name][atom_index])
    with open('coord.crdchg', 'w') as f:
        f.write(snap_shot)

def from_xyz_to_crdchg(nmol, natoms, mol_list, dic_mol_natom,dic_mol_element,dic_mol_chg,filename='coord.xyz'):
    with open(filename, 'r') as f:
        lines = f.readlines()
    current_mol = 0
    mol_num = []
    
    mol_num = [0]
    for _,i in enumerate(mol_list):
        if _ == 0:
            mol_num.append(dic_mol_natom[i])
        else:
            mol_num.append(dic_mol_natom[i]+mol_num[_])
    tot_natom = mol_num[-1]
    tot_mol = int(len(lines)/(tot_natom+2))
    #print(tot_mol)
    resnum_old = 0
    for current_t_mol in range(tot_mol):
        snap_shot = "{nmol}  {natoms}\n-----\n".format(nmol=nmol,natoms=natoms)
        current_mol = 0
        for index,line in enumerate(lines[current_t_mol*(natoms+2)+2:(current_t_mol+1)*(natoms+2)]):
            while True:
                if index < mol_num[current_mol+1]: # this is the correct molecular.
                    atom_index = index - mol_num[current_mol]
                    break
                else:
                    current_mol += 1
            mol_name = mol_list[current_mol]
            if atom_index == 0:
                if mol_name=='Wat':snap_shot += '0 1 {index}\n'.format(index=dic_mol_index[mol_name])
                elif mol_name=='Azo':snap_shot += '0 1 {index}\n'.format(index=dic_mol_index[mol_name])
                else: snap_shot += '0 1 {index}\n'.format(index=dic_mol_index[mol_name])
                
            snap_shot += "{ele:4s}{x:15.8f}{y:15.8f}{z:15.8f}{chg:11.5f}\n".format(ele=dic_mol_element[mol_name][atom_index],
                                                                               x=float(line.split()[1]),
                                                                               y=float(line.split()[2]),
                                                                               z=float(line.split()[3]),  
                                                                               chg=dic_mol_chg[mol_name][atom_index])
        with open('cluster_{}.crdchg'.format(current_t_mol), 'w') as f:
            f.write(snap_shot)


if __name__ == "__main__":
    nmol = 6
    natoms = 25
    Azo_chg = [-0.101115,-0.105311,-0.372597, 0.156724, 0.161005, 0.156522, -0.369976, 0.156872, 0.160781, 0.157095]
    wat_chg = [-0.834, 0.417, 0.417]
    mol_list= ['Azo',
               'Wat','Wat','Wat','Wat','Wat']
              # 'Wat','Wat','Wat','Wat','Wat',
              # 'Wat','Wat','Wat','Wat']
    dic_mol_natom = {'Azo':10,'Wat':3}
    dic_mol_index = {'Azo':1, 'Wat':2}
    dic_mol_element={'Azo':['N','N','C','H','H','H','C','H','H','H'],
                     'Wat':['O','H','H']}
    dic_mol_chg   = {'Azo':Azo_chg, 'Wat':wat_chg}
    filename = '../../GTSH.XYZ'
    #print(filename.split('.')[-1])
    if filename.split('.')[-1] == 'pdb':
        from_pdb_to_crdchg(filename,nmol,natoms,Azo_chg,wat_chg)
    elif filename.split('.')[-1] == 'txt':
        from_txt_to_crdchg(nmol,natoms,mol_list,
                           dic_mol_natom,dic_mol_element,dic_mol_chg)
    elif filename.split('.')[-1].lower() == 'xyz':
        from_xyz_to_crdchg(nmol,natoms,mol_list,
                           dic_mol_natom,dic_mol_element,dic_mol_chg,filename=filename)
