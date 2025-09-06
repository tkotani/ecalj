#!/usr/bin/env python3
import sys, os, argparse
'''
This script refines a crystal structure from a POSCAR file using symmetry analysis.
It reads a POSCAR file, analyzes its symmetry, and outputs a refined structure to a new POSCAR file.
The script uses pymatgen's SpacegroupAnalyzer to refine the structure based on symmetry operations.
'''

### Options
parser = argparse.ArgumentParser(description="Refine crystal structure from POSCAR file")
parser.add_argument('--poscar', type=str, required=True, help='Input POSCAR file path')
parser.add_argument('--refine', type=str, required=True, help='Output refined POSCAR file path')
parser.add_argument('--symprec1', type=float, default=3e-3, help='Initial symmetry precision (default: 1e-5)')
parser.add_argument('--symprec2', type=float, default=1e-8, help='Final symmetry precision (default: 1e-6)')
parser.add_argument('--primitive', action='store_true', help='Output primitive standard structure instead of refined structure')
parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
args = parser.parse_args()

def refine_structure(poscar_file, output_file, symprec1=3e-3, symprec2=1e-8, primitive=False, verbose=False):
    """
    Refine crystal structure from POSCAR file
    
    Parameters:
    poscar_file: Input POSCAR file path
    output_file: Output refined POSCAR file path
    symprec1: Initial symmetry precision for analysis
    symprec2: Final symmetry precision for refinement
    primitive: If True, output primitive standard structure
    verbose: If True, print detailed information
    """
    
    from pymatgen.core import Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    
    try:
        # Load POSCAR file
        if verbose:
            print(f"Loading structure from: {poscar_file}")
        struc = Structure.from_file(poscar_file)
        
        # Analyze initial structure
        if verbose:
            print(f"Initial structure analysis with symprec={symprec1}")
        init_analyzer = SpacegroupAnalyzer(struc, symprec=symprec1)
        init_sg = init_analyzer.get_space_group_number()
        init_symbol = init_analyzer.get_space_group_symbol()
        
        if verbose:
            print(f"Initial space group: {init_sg} ({init_symbol})")
            print(f"Initial lattice parameters: a={struc.lattice.a:.6f}, b={struc.lattice.b:.6f}, c={struc.lattice.c:.6f}")
            print(f"Initial angles: α={struc.lattice.alpha:.2f}°, β={struc.lattice.beta:.2f}°, γ={struc.lattice.gamma:.2f}°")
            print(f"Initial number of sites: {len(struc.sites)}")
        
        # Refine structure with higher precision
        if verbose:
            print(f"Refining structure with symprec={symprec2}")
        refined_struc = init_analyzer.get_refined_structure()
        final_analyzer = SpacegroupAnalyzer(refined_struc, symprec=symprec2)
        
        # Choose output structure type
        if primitive:
            output_struc = final_analyzer.get_primitive_standard_structure()
            if verbose:
                print("Using primitive standard structure")
        else:
            output_struc = final_analyzer.get_refined_structure()
            if verbose:
                print("Using refined structure")
        
        # Verify final space group
        final_sg = final_analyzer.get_space_group_number()
        final_symbol = final_analyzer.get_space_group_symbol()
        
        if verbose:
            print(f"Final space group: {final_sg} ({final_symbol})")
            print(f"Final lattice parameters: a={output_struc.lattice.a:.6f}, b={output_struc.lattice.b:.6f}, c={output_struc.lattice.c:.6f}")
            print(f"Final angles: α={output_struc.lattice.alpha:.2f}°, β={output_struc.lattice.beta:.2f}°, γ={output_struc.lattice.gamma:.2f}°")
            print(f"Final number of sites: {len(output_struc.sites)}")
        
        # Save refined structure to output file
        output_struc.to(fmt="poscar", filename=output_file)
        print(f"Refined structure saved to: {output_file}")
        
        # Check if space group changed
        if init_sg != final_sg:
            print(f"Warning: Space group changed from {init_sg} ({init_symbol}) to {final_sg} ({final_symbol})")
        else:
            print(f"Space group maintained: {final_sg} ({final_symbol})")
            
        # Print composition
        print(f"Composition: {output_struc.composition.reduced_formula}")
        
        return True
        
    except FileNotFoundError:
        print(f"Error: Input file '{poscar_file}' not found.")
        return False
    except Exception as e:
        print(f"Error processing structure: {str(e)}")
        return False

def main():
    """Main function"""
    
    # Check if input file exists
    if not os.path.exists(args.poscar):
        print(f"Error: Input file '{args.poscar}' does not exist.")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.refine)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        if args.verbose:
            print(f"Created output directory: {output_dir}")
    
    # Refine structure
    success = refine_structure(
        poscar_file=args.poscar,
        output_file=args.refine,
        symprec1=args.symprec1,
        symprec2=args.symprec2,
        primitive=args.primitive,
        verbose=args.verbose
    )
    
    if success:
        print("Structure refinement completed successfully.")
        sys.exit(0)
    else:
        print("Structure refinement failed.")
        sys.exit(1)

if __name__ == "__main__":
    main()
