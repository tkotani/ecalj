#!/usr/bin/env python3
import sys, os, configparser, argparse
from datetime import datetime
'''
This script queries the Materials Project database for materials with specific properties and saves the results to a file.
It allows users to specify various criteria such as the number of sites, elements, and other properties.
The script uses the MPRester class from the mp_api.client module to interact with the Materials Project API.
It also includes options for filtering materials based on their magnetic properties and whether they are theoretical or experimental.
The script can save the queried materials in POSCAR format and includes options for refining the structure based on symmetry.'''


### Options
parser = argparse.ArgumentParser(description="Get Targets from Materials Project")
parser.add_argument('--fname', type=str, default="joblist")
parser.add_argument('--nsites', nargs='+', type=int, default=[1, 6], help='Number of sites or range (min, max)')
parser.add_argument('--nels', nargs='+', type=int, default=None)
parser.add_argument('--include', nargs='+', type=str, default=['ANY'], help='Required elements like Si O .., except excluded atoms. ANY(default) is to all')
parser.add_argument('--exclude', nargs='+', type=str, default=['NG','LN','AC'])
parser.add_argument('--mag', type=str, choices=['true','false','both'], default='false')
parser.add_argument('--theoretical', type=str, choices=['true','false','both'], default='false')
parser.add_argument('--metal', type=str, choices=['true','false','both'], default='false')
parser.add_argument('--mpid', type=str, nargs='+', default=False, help='Enter material_ids like 1 10 100 .... To obtain a dataset for comparison with experimental values, enter "expt"')
parser.add_argument('--poscar', type=bool, default=True)
parser.add_argument('--dir', type=str, default='lists')
parser.add_argument('--conf', type=str, default='./config.ini')
args = parser.parse_args(sys.argv[1:])

print('loading MPRester ...')
from mp_api.client import MPRester
### Get API-KEY from configuration file
config = configparser.ConfigParser()
config.read(args.conf)
apikey = config['DEFAULT'].get('apikey')

### Set elements: NobleGas, Lanthanoide, Actinoide
special_elements = {
    'NG' : 'He Ne Ar Xe Kr Xe Rn',
    'LN' : 'La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu',
    'AC' : 'Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg'
}
expt_id = '10044 10695 1070 1078292 1087 1132 1156 11674 1190 1201 1239 12779 1317 1321 1367 14 14091 14092 1479 149 1541 1550 15700 15776 15777 1836 19 19006 19009 1943 1944 19717 19765 19833 20012 20194 20305 20351 20485 20554 20724 20731 20832 21276 2133 21374 2172 2176 2201 22205 22261 22304 22607 22736 22811 22894 22895 22913 22914 22925 2343 2395 2490 252 2534 2624 3078 32 3345 34202 3497 3595 3668 3772 3829 3839 4008 404 406 408 4175 422 4452 452 4524 4666 4730 4756 4763 4809 4840 4899 4953 4979 5190 5213 5342 5350 541837 542812 5518 556911 5702 672 675626 704645 7475 753 7631 8016 8017 804 8062 8182 856 9814'

try:
    if args.mpid[0] == "expt":
        mpid = expt_id.split()
    else:
        mpid = args.mpid
except: pass


class Conditions:

    def __init__(self):
        
        self.query_dict = {
            "fields": [
                "material_id",
                "nsites",
                "elements",
                "band_gap",
                "symmetry",
                "formula_pretty",
                "structure",
            ]
        }
        if args.mpid:
            mpid_list = [
                f"mp-{item}" if item.isdigit() else item
                for item in mpid
                if item.isdigit() or item.startswith("mp-")
            ]
            self.query_dict["material_ids"] = mpid_list
        self.n_sites()
        self.n_els()
        self.get_flags()
        self.required = self.target_elements(args.include)
        self.avoided = self.target_elements(args.exclude)
        
    def n_sites(self):

        if len(args.nsites) == 1:
            self.nsites = tuple(args.nsites * 2)
        else:
            self.nsites = tuple(args.nsites)
        self.query_dict["num_sites"] = self.nsites

    def n_els(self):

        if args.nels is not None:
            if len(args.nels) == 1:
                self.nels = tuple(args.nels * 2)
            else:
                self.nels = tuple(args.nels)
            self.query_dict["num_elements"] = self.nels

    def get_flags(self):

        [mag, theoretical, metal] = [self.str_to_bool(getattr(args, v)) for v in ['mag', 'theoretical', 'metal']]
        for flag,v in zip([theoretical, metal], ['theoretical', 'is_metal']):
            if flag is None: continue
            self.query_dict[v] = flag

        if mag is not None:
            self.query_dict["fields"].append("is_magnetic")
        #if mag == False:
            #self.query_dict["total_magnetization"] = (0, 0.01)
        self.mag = mag

    @staticmethod
    def str_to_bool(value):

        if value.lower() in ('true', '1'):
            return True
        elif value.lower() in ('false', '0'):
            return False
        elif value.lower() == 'both':
            return None
        else:
            raise argparse.ArgumentTypeError(f"Boolean value expected, got {value}")


    @staticmethod
    def target_elements(condition):

        for value in condition:
            if value.lower() in ('any', 'none'):
                return value.lower()
        
        els = ''
        for el in condition:
            els += ' ' + special_elements[el] if el in ['NG', 'LN', 'AC'] else ' ' + el
        elements = els.split()
        return elements


def get_magnetic(struc):

    from pymatgen.analysis.magnetism.analyzer import CollinearMagneticStructureAnalyzer as CMSA
    try:
        st = CMSA(struc)
        msm = st.magnetic_species_and_magmoms
        ordering = str(st.ordering).split('.')[1]
        mag_info = f"{msm} {ordering}"
    except:
        mag_info = "unknown"
    return mag_info


def get_list(dir_name):

    ### get conditions
    conditions = Conditions()
    query_dict = conditions.query_dict
    num_sites = conditions.nsites
    required = conditions.required
    avoided = conditions.avoided
    mag = conditions.mag

    with open(f'{dir_name}/datearg.chk', 'w') as f:
        xdate = datetime.now().strftime("%Y%m%d-%H%M%S")
        print(f'# date: {xdate}',file=f)
        print(f'# args: {args}',file=f)

    with MPRester(apikey) as mpr:

        docs=mpr.materials.summary.search( **query_dict )
        if args.poscar: os.makedirs(f'{dir_name}/POSCARALLnoREFINE',exist_ok=True)
        if args.poscar: os.makedirs(f'{dir_name}/POSCARALL',exist_ok=True)
        ff=open(f'{dir_name}/{args.fname}','w')
        print('# material_id nsites formula bandgap symmetry-number symmetry-symbol magnetic', file=ff, end='\n')
        count=0
        for doc in docs:

            ### Check if elements are consistent with Targets or not
            if (
                    mag not in [doc.is_magnetic, None] or
                    (avoided != "none" and any(str(element) in avoided for element in doc.elements)) or
                    (required != "any" and not any(str(element) in required for element in doc.elements))
            ):
                continue

            # save to poscar
            if args.poscar:
                doc.structure.to(fmt='poscar', filename=f'{dir_name}/POSCARALLnoREFINE/POSCAR.{doc.material_id}')
                refine_structure(doc.structure,doc.symmetry.number, out=f'{dir_name}/POSCARALL/POSCAR.{doc.material_id}')
            if doc.is_magnetic == False:
                mag_info = "NM"
            else:
                mag_info = get_magnetic(doc.structure)

#            print(doc.material_id.strip('mp-'),
            print(doc.material_id,
                  doc.nsites,
                  doc.formula_pretty,
                  '{:.04f}'.format(doc.band_gap),
                  doc.symmetry.number,
                  doc.symmetry.symbol,
                  mag_info,
                  file=ff)
            count+=1

    ff.close()
    print('#####   file       :', dir_name + '/' + args.fname)
    print('  number of sites  :', num_sites)
    print('  number of mp     :', count)
    
from pymatgen.core import Structure
from pymatgen.symmetry.groups import SpaceGroup
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
def refine_structure(struc, sgn, out, symprec1=1e-5, symprec2=1e-6):
    # Load POSCAR & space group
    #struc = Structure.from_file(poscar)
    sg = SpaceGroup.from_int_number(sgn)

    # Analyze initial structure (precision=1e-3)
    symprec0 = symprec1
    while True:
        init_analyzer = SpacegroupAnalyzer(struc, symprec=symprec0)
        init_sg = init_analyzer.get_space_group_number()
        if init_sg == sgn:
            break
        print('Initial detected space group:', init_sg)
        symprec0 *= 10

    # Refine structure (precision=1e-6)
    refined_struc = init_analyzer.get_refined_structure()
    #print(refined_struc)
    final_analyzer = SpacegroupAnalyzer(refined_struc, symprec=symprec2)
    corrected_struc = final_analyzer.get_refined_structure()

    # Save structure to POSCAR
    primitive_struc=final_analyzer.get_primitive_standard_structure()
    primitive_struc.to(fmt="poscar", filename=out)
    print(primitive_struc)
    #corrected_struc.to(fmt="poscar", filename=out)
    print(f"Refined structure saved to {out}")

    # Verify the space group of the refined structure
    final_sg = final_analyzer.get_space_group_number()
    print(f"Final detected space group: {final_sg}")

    if final_sg != sgn: #space_group_number:
        print("Warning: Refined structure does not match the specified space group.")
    else:
        print("Refined structure matches the specified space group.")

if __name__=="__main__":
    print('Start query to mp')
    dir_name = f"INPUT/{args.dir}"
    os.makedirs(dir_name, exist_ok=True)
    get_list(dir_name)
