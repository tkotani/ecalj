from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.magnetism.analyzer import CollinearMagneticStructureAnalyzer as CMSA
import re, shutil

element_data = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20,
    'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27,'Ni': 28, 'Cu': 29, 'Zn': 30,
    'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40,
    'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
    'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60,
    'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70,
    'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80,
    'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr':87, 'Ra': 88, 'Ac': 89, 'Th': 90,
    'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100,
    'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110,
    'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118
}

s_block = ['H', 'He', 'Li', 'Be', 'Na', 'Mg', 'K', 'Ca', 'Rb', 'Sr', 'Cs', 'Ba', 'Fr', 'Ra']
p_block = ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']
d_block = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Lu', 'Lr']
f_block = ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']

class SetMag:
    '''   Set the magnetic moment of each atom in the POSCAR file.
    The magnetic moment is set in the following order:
    1. Ferromagnetic (FM)
    2. Ferrimagnetic (FiM)
    3. Antiferromagnetic (AFM)
    4. Non-magnetic (NM)'''
    def __init__(self, key, num, poscar, nr, ordering):

        self.num = num
        self.nr = nr
        self.ordering = ordering

        with MPRester(key) as mpr:
            struc = mpr.get_structure_by_material_id(num)
        st = CMSA(struc)
        self.msm = st.magnetic_species_and_magmoms
        if ordering == "AFM":
            self.magmoms = struc.site_properties['magmom']

        with open(poscar,'r') as pos:
            elements = re.sub(r'[0-9]+','', pos.readline())
            el_list = elements.split()
        self.el_z = { el: element_data[el] for el in el_list}

        self.get_mag()

        
    def get_mag(self):

        if self.ordering in ["FM", "FiM"]:
            for atom, mmom in self.msm.items():
                if isinstance(mmom, list):
                    self.msm[atom] = sum(mmom) / len(mmom)
            self.fctrl = open(f"ctrls.{self.num}", "a")
            self.fctrl.writelines('SPEC\n')

        if self.ordering == "FM":
            self.process_ferromagnetic()
        elif self.ordering == "FiM":
            self.process_ferrimagnetic()
        elif self.ordering == "AFM":
            self.process_antiferromagnetic()

        for atom, z in self.el_z.items():
            self.print_mmom(atom, z, 0)
        self.fctrl.close()


    def process_ferromagnetic(self):

        for atom, mmom in self.msm.items():
            z = self.el_z[atom]
            self.print_mmom(atom, z, mmom)
            del self.el_z[atom]


    def process_ferrimagnetic(self):

        mags = [mmom for mmom in self.msm.values() if mmom >= 0.1]
        for atom, mmom in self.msm():
            if len(mags) == 0: break
            if len(mags) > 1 and mmom == min(mags):
                mmom = -mmom
            z = self.el_z[atom]
            self.print_mmom(atom, z, mmom)
            del self.el_z[atom]

    def process_antiferromagnetic(self):

        shutil.copy(f'ctrls.{self.num}', 'ctrls0')
        self.fctrl = open(f"ctrls.{self.num}", 'w')
        with open(f"ctrls0", "r") as original:
            i = 0
            for line in original:
                if 'STRUC' in line or 'ATOM=' not in line:
                    self.fctrl.write(line)
                    continue
                elif 'ATOM=' in line:
                    atom = (line.split()[0]).split('ATOM=')[1]
                    if abs(self.magmoms[i]) < 0.01:
                        self.fctrl.write(line)
                    elif self.magmoms[i] > 0:
                        self.fctrl.write(re.sub(atom, atom + 'up', line))
                    else:
                        self.fctrl.write(re.sub(atom, atom + 'dn', line))
                    i += 1

        self.fctrl.writelines('SPEC\n')
        for atom, magmom in zip(self.msm.keys(), self.magmoms):
            if abs(magmom) < 0.01:
                continue
            elif magmom > 0:
                z = self.el_z[atom]
                self.print_mmom(f"{atom}up", z, magmom)
            else:
                z = self.element_z[atom]
                self.print_mmom(f"{atom}dn", z, magmom)
            del self.el_z[atom]


    def print_mmom(self, atom, z, mmom):

        if atom in s_block: mmm = f"{mmom} 0 0 0"
        elif atom in p_block: mmm = f"0 {mmom} 0 0"
        elif atom in d_block: mmm = f"0 0 {mmom} 0"
        elif atom in f_block: mmm = f"0 0  0 {mmom}"
        spec = f"     ATOM={atom} Z={z} MMOM={mmm}"
        print(spec, end='\n', file=self.fctrl)
