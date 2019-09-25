import numpy as np
from dpdata.periodic_table import Element

def _make_siesta_01_common(sys_data, fp_params):
    tot_natoms = sum(sys_data['atom_numbs'])
    ntypes = len(sys_data['atom_names'])

    ret = ""
    ret += 'SystemName        system\n'
    ret += 'SystemLabel       system\n'
    ret += 'NumberOfAtoms     %d\n' % tot_natoms
    ret += 'NumberOfSpecies   %d\n' % ntypes
    ret += '\n'
    ret += 'WriteForces       T\n'
    ret += 'WriteCoorStep     T\n'
    ret += 'WriteCoorXmol     T\n'
    ret += 'WriteMDXmol       T\n'
    ret += 'WriteMDHistory    T\n\n'

    if 'ecut' in fp_params.keys():
        ecut = fp_params['ecut']
        ret += 'MeshCutoff            %s' % str(ecut)
        ret += ' Ry\n'
    if 'ediff' in fp_params.keys():
        ediff = fp_params['ediff']
        ret += 'DM.Tolerance          %e\n' % ediff
    if 'mixWeight' in fp_params.keys():
        mixingWeight = fp_params['mixingWeight']
        ret += 'DM.MixingWeight       %f\n' % mixingWeight
    if 'NumberPulay' in fp_params.keys():
        NumberPulay = fp_params['NumberPulay']
        ret += 'DM.NumberPulay         %d\n' % NumberPulay
    ret += 'DM.UseSaveDM          true\n'
    ret += 'XC.functional          GGA\n'
    ret += 'XC.authors             PBE\n'
    ret += 'MD.UseSaveXV           T\n\n'
    ret += 'DM.UseSaveDM           F\n'
    ret += 'WriteDM                F\n'
    ret += 'WriteDM.NetCDF         F\n'
    ret += 'WriteDMHS.NetCDF       F\n'
    return ret

def _make_siesta_02_species(sys_data, pps):
    atom_nums = sys_data['atom_numbs']
    atom_names = sys_data['atom_names']
    ntypes = len(atom_nums)
    assert (ntypes == len(atom_names))
    assert (ntypes == len(pps))
    ret = ''
    ret += '%block Chemical_Species_label\n'
    for i in range(0, len(atom_names)):
        ret += str(i + 1) + '\t' + str(Element(atom_names[i]).Z) + '\t' + atom_names[i] + '\n'
    ret += '%endblock Chemical_Species_label\n'
    return ret

# ## kpoints !!!
def _make_siesta_03_kpoint(sys_data, fp_param):
    if 'kspacing' in fp_param.keys():
        kspacing = fp_param['kspacing']
        cell = sys_data['cells'][0]
        cell = np.reshape(cell, [3, 3])

        rcell = np.linalg.inv(cell)

        rcell = rcell.T

        kpoints = [(np.ceil(2 * np.pi * np.linalg.norm(ii) / kspacing).astype(int))
                   for ii in rcell]

        ret = ""
        ret += '%block kgrid_Monkhorst_Pack\n'
        ret += '%d' % kpoints[0]
        ret += '\t0\t0\t0.0\n'

        ret += '0\t'
        ret += '%d' % kpoints[1]
        ret += '\t0\t0.0\n'

        ret += '0\t0\t'
        ret += '%d' % kpoints[2]
        ret += '\t0.0\n'

        ret += '%endblock kgrid_Monkhorst_Pack\n'
        return ret
    else:
        return ''

### coordinate
def _make_siesta_04_ucVectorCoord(sys_data):
    cell = sys_data['cells'][0]
    cell = np.reshape(cell, [3, 3])
    coordinates = sys_data['coords'][0]
    atom_names = (sys_data['atom_names'])
    atom_numbs = (sys_data['atom_numbs'])
    ntypes = len(atom_names)
    ret = ""
    ret += "LatticeConstant 1.00 Ang\n"
    ret += "%block LatticeVectors\n"
    for ii in range(3):
        for jj in range(3):
            ret += "%f " % cell[ii][jj]
        ret += "\n"
    ret += "%endblock LatticeVectors\n"

    ret += "\n"
    ret += "AtomicCoordinatesFormat Ang\n"
    ret += "%block AtomicCoordinatesAndAtomicSpecies\n"
    cc = 0
    for ii in range(ntypes):
        for jj in range(atom_numbs[ii]):
            ret += "%f %f %f %d %s\n" % (coordinates[cc][0],
                                         coordinates[cc][1],
                                         coordinates[cc][2],
                                         ii + 1,
                                         atom_names[ii])
            cc += 1
    ret += "%endblock AtomicCoordinatesAndAtomicSpecies"
    return ret

def make_siesta_input(sys_data, fp_pp_files, fp_params):
    ret = ""
    ret += _make_siesta_01_common(sys_data, fp_params)
    ret += "\n"
    ret += _make_siesta_02_species(sys_data, fp_pp_files)
    ret += "\n"
    ret += _make_siesta_03_kpoint(sys_data, fp_params)
    ret += "\n"
    ret += _make_siesta_04_ucVectorCoord(sys_data)
    ret += "\n"
    return ret

# sys_data = {'atom_numbs':[4, 1], 'atom_names':['H','C'], 'cells':[[10,0,0,0,10,0,0,0,10]], \
#             'mixingWeight': 0.05, 'NumberPulay': 5,
#             'coords':[[[5.36586,4.05415,3.59848], [3.92611,4.82628,4.42448], [5.5636,5.65033,4.43064], \
#                      [5.28729,4.14347,5.36703],[5.01936,4.68879,4.44103]]]}
# jdata = {'fp_pp_files':["C.psf", "H.psf"], 'fp_params':{'ecut': 300, 'ediff': 1.0e-4, 'kspacing': 2, \
#          'mixingWeight': 0.05, 'NumberPulay':5 }}
# fp_pp_files = jdata['fp_pp_files']
# fp_params = jdata['fp_params']
# ret = make_siesta_input(sys_data,fp_pp_files, fp_params)
# f = open('input', 'w')
# f.write(ret)
# f.close()

#############################read output#####################################
# def get_single_line_tail(fin, keyword):
#     file = open(fin, 'r')
#     res = []
#     for value in file:
#         if keyword in value:
#             temp = len(value.split()) - 1
#             res.append(float(value.split()[temp]))
#             file.close()
#             return res
#     return res

## atomnum: number of atoms,  row numbers
## begin_column: begin column num
## column_num: read column num
# def extract_keyword(fout, keyword, down_line_num, begin_column, column_num):
#     file = open(fout, 'r')
#     ret = []
#     flag = 0
#     idx = 0
#     # for (num,value) in enumerate(file):
#     for value in file:
#         if keyword in value:
#             flag = 1
#             continue
#         if flag == 1:
#             if idx < down_line_num:
#                 idx += 1
#             else:
#                 flag = 0
#                 continue
#             for i in range(begin_column, column_num):
#                 if not value.split()[i].isalpha():
#                     ret.append(float(value.strip().split()[i]))
#                 else:
#                     ret.append(value.strip().split()[i])
#             continue
#     file.close()
#     return ret
#
# def get_atom_types(fout, atomnums):
#     covert_type = extract_keyword(fout, 'outcoor: Atomic coordinates (Ang):', atomnums, 3, 4)
#     atomtype = []
#     for i in range(0, len(covert_type)):
#         atomtype.append(int(covert_type[i]) - 1)
#     return atomtype
#
# def get_virial(fout, cells):
#     vols = []
#     for ii in cells:
#         ### calucate vol
#         vols.append(np.linalg.det(ii.reshape([3, 3])))
#     ret = extract_keyword(fout, 'siesta: Stress tensor (static) (eV/Ang**3):', 3, 1, 4)
#     ret = np.array([ret])
#     for idx, ii in enumerate(ret):
#         ## siesta: 1eV/A^3= 1.60217*10^11 Pa ,  ---> qe: kBar=10^6Pa  --> 1 Bar = 10^5Pa	1KBar = 10^8 Pa
#         ii *= vols[idx] * 1e3 / 1.602176621e6 * (1.602176621e3)
#     return ret
#
# def get_atom_numbs(atomtypes):
#     atom_numbs = []
#     for i in set(atomtypes):
#         atom_numbs.append(atomtypes.count(i))
#     return atom_numbs

# def get_frame(fin, fout):
#     tot_natoms = int(get_single_line_tail(fin, 'NumberOfAtoms')[0])
#     NumberOfSpecies = int(get_single_line_tail(fin, 'NumberOfSpecies')[0])
#     atom_names = extract_keyword(fin, '%block Chemical_Species_label', NumberOfSpecies, 2, 3)
#
#     atom_types = get_atom_types(fout, tot_natoms)
#     atom_numbs = get_atom_numbs(atom_types)
#
#     cells = extract_keyword(fout, 'outcell: Unit cell vectors (Ang):', 3, 0, 3)
#     coords = extract_keyword(fout, 'outcoor: Atomic coordinates (Ang):', tot_natoms, 0, 3)
#
#     energies = get_single_line_tail(fout, 'siesta: E_KS(eV) =')
#     forces = extract_keyword(fout, 'siesta: Atomic forces (eV/Ang):', tot_natoms, 1, 4)
#
#     data = {}
#     data['orig'] = np.array([0, 0, 0])
#     data['atom_names'] = atom_names
#     data['atom_numbs'] = atom_numbs
#     data['atom_types'] = np.array(atom_types)
#     data['coords'] = np.array([coords])
#     data['cells'] = np.array([cells])
#     data['energies'] = np.array([energies])
#     data['forces'] = np.array([forces])
#     data['virials'] = get_virial(fout, data['cells'])
#     return data

# data = get_frame('input', 'output')
