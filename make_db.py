import pandas as pd
import pymol
from pymol import cmd
import os
import argparse
import numpy
from openmm.app import PDBFile
from pdbfixer.pdbfixer import PDBFixer, proteinResidues, dnaResidues, rnaResidues, _guessFileFormat
import openmm as mm
import MDAnalysis as mda

parser = argparse.ArgumentParser()

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("input", type=str, help="directory with .pdb or .cif files OR text file containing PDB codes (one per line)")
parser.add_argument("outputDir", type=str, help="name of the directory to save results to (structures and dataset)")
parser.add_argument("-t", "--threshold", default="15", type=int, help="maximum length (in residues) of biological polymer to keep, default = 15")
parser.add_argument("-c", "--clean", action="store_true", help="if true, clean up files by removing extra chains")
parser.add_argument("-s", "--split", action="store_true", help="if true, split receptor and peptide")
parser.add_argument("-b", "--bsa", action="store_true", help="if true, calculate SASA and BSA")
parser.add_argument("-p", "--peptide", action="store_true", help="if true, only keep biological polymers - delete small molecules, solvents, and ions")
args = parser.parse_args()

print(args.threshold)

dirName = args.outputDir

if not os.path.exists(dirName):
    os.mkdir(dirName)

#######################################################################################################################################

def get_pdb_files():
    
    df = pd.read_csv(args.input, sep=',',header=None)
    pdb_codes = df.to_numpy()
    os.chdir(dirName)
    os.mkdir("base_pdbs")
    os.mkdir("isolated_pdbs")
    os.mkdir("fixed_isolated_pdbs")
    os.mkdir("solvated_fixed_isolated_pdbs")
    os.mkdir("receptors")
    os.mkdir("peptides")
    os.mkdir("temp")
    os.chdir("temp")

    for entry in pdb_codes:
        try:
            entry = (', '.join(entry))
            print("Fetching PDB ID", entry)
            cmd.reinitialize()
            cmd.fetch(entry)
            cmd.remove("not (alt ""+A)")
            cmd.alter('all', "alt=''")
            cmd.save(os.path.join(path, "base_pdbs", entry+".pdb"))
        except:
            print("failed at getting PDB ID", entry)
            
    
    os.chdir("..")
    os.chdir("..")

###########################################################################################################################################

def takeSecond(lst):
    return lst[1]

def pdb_separator():

    # loop through files and calculates SASAs and BSA #
    for filename in os.listdir(os.path.join(path, "base_pdbs")):
        if filename.endswith(".pdb") or filename.endswith(".cif"):
            filenamefull = os.path.join(path, "base_pdbs", filename)
            chain_lens = []
            chain_lens_ca = []
            i = 0
            repetition_check = 0
        #try:
            # reloads pymol and loads file #
            cmd.reinitialize()
            cmd.load(filenamefull)

            # gets PDB ID without extension #
            file = os.path.splitext(filename)

            # gets all chains #
            chains = cmd.get_chains(file[0])

            print("Cleaning chains for file ", file[0])

            # looks for the chain with most atoms (assumes this is the most complete chain of the protein) #
            # then, delete all other chains that have less than X alpha carbons (thus, keeps peptides)    #
            for chain in chains:
                count = cmd.count_atoms(file[0] + ' and chain ' + chain)
                count_ca = cmd.count_atoms(file[0] + ' and chain ' + chain + ' and name CA')
                chain_lens.append(count)
                chain_lens_ca.append(count_ca)

            for n in chain_lens:
                if n < max(chain_lens) and chain_lens_ca[i] > args.threshold:
                    cmd.remove("chain "+chains[i])
                    i = i + 1
                elif n == max(chain_lens) and repetition_check == 0:
                    repetition_check = 1
                    i = i + 1
                elif n == max(chain_lens) and repetition_check == 1:
                    cmd.remove("chain "+chains[i])
                    i = i + 1
                else:
                    i = i + 1

            max_value = max(chain_lens)
            max_index = chain_lens.index(max_value)

            # gets the peptide that is interacting with the leftover chain, selects the complex and deletes the rest #
            if args.peptide:
                cmd.create("complex", (file[0] + ' and chain '+chains[max_index]+' + bychain all within 3 of chain '+chains[max_index]))
                cmd.create("complex", "complex and polymer")
            else:
                cmd.create("complex", (file[0] + ' and chain '+chains[max_index]+' + bymolecule all within 3 of chain '+chains[max_index]))
            cmd.remove("not complex")
            filename = 'clean_' + file[0]+".pdb"
			
			
			# reorders chains so the receptor comes first as chain A and the peptide follows it as chain B
            chains = cmd.get_chains("complex")

            lst = []
            for chain in chains:
                count_ca = cmd.count_atoms('complex and chain ' + chain + ' and name CA')
                lst.append([chain, count_ca])
            
            lst.sort(reverse=False,key=takeSecond)

            placeholder = ['ZZ', 'YY', 'XX', 'VV', 'WW', 'FF', 'TT', 'PP', 'QQ']
            i = 0
            for chain in chains:
                cmd.alter("chain "+ chain , "chain='"+placeholder[i]+"'" )
                #cmd.set("retain_order",0)
                i = i+1
            
            i = 0
            chains = cmd.get_chains("complex")
            for chain in chains:
                cmd.alter("chain "+ chain , "chain='"+lst[i][0]+"'" )
                #cmd.set("retain_order",1)
                i = i+1
            
            cmd.save(os.path.join(path, "isolated_pdbs", filename),"complex",-1)

        #except:
            #print("failed at cleaning PDB ID", filename) 

###########################################################################################################################################  

def process_files():
        
	# fills in missing residues, atoms, adds hydrogens, chooses the highest occupancy alternative residue position (if theres one), and solvates
    for filename in os.listdir(os.path.join(path, "isolated_pdbs")):
        file = os.path.splitext(filename)
        filenamefull = os.path.join(path, "isolated_pdbs", filename)
        print("Processing file ", file[0])
        fixer = PDBFixer(filenamefull)
        fixer.findMissingResidues()
        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        for key in keys:
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)

        PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(path, "fixed_isolated_pdbs", file[0]+"_processed.pdb"), 'w'))

        fixer.removeHeterogens(True)
        maxSize = max(max((pos[i] for pos in fixer.positions))-min((pos[i] for pos in fixer.positions)) for i in range(3))
        boxSize = maxSize*mm.Vec3(1, 1, 1)
        print("Solvating file ", file[0])
        fixer.addSolvent(boxSize)

        PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(path, "solvated_fixed_isolated_pdbs", file[0]+"_processed.pdb"), 'w'))
    

#################################################################################################################################################	

def split_files():
    
	# for each complex, create split receptor/peptide files and save them along with their .fasta sequences
    for filename in os.listdir(os.path.join(path, "solvated_fixed_isolated_pdbs")):
        lst = []
        if filename.endswith(".pdb") or filename.endswith(".cif"):
            filenamefull = os.path.join(path, "solvated_fixed_isolated_pdbs", filename)
            cmd.reinitialize()
            cmd.load(filenamefull)
            file = os.path.splitext(filename)
            # gets first two chains after removing duplicate entries #
            chains = cmd.get_chains(file[0])

            for chain in chains:
                count_ca = cmd.count_atoms(file[0]+' and chain ' + chain + ' and name CA')
                lst.append([chain, count_ca])

            lst.sort(reverse=False,key=takeSecond)
            placeholder = ['ZZ', 'YY', 'XX', 'VV', 'WW', 'FF', 'TT', 'PP', 'QQ']
            
            i = 0
            for chain in chains:
                cmd.alter("chain "+ chain , "chain='"+placeholder[i]+"'" )
                i = i+1

            i = 0
            chains = cmd.get_chains(file[0])
            for chain in chains:
                cmd.alter("chain "+ chain , "chain='"+lst[i][0]+"'" )
                i = i+1
            
            
            chains = cmd.get_chains(file[0])
            # create objects for protein, peptide, and protein+peptide pair #
            cmd.create("r_"+file[0], (file[0] + ' and chain '+chains[0]))
            cmd.create("p_"+file[0], (file[0] + ' and chain '+chains[1]))
            
            f_protein = str(cmd.get_fastastr("r_"+file[0]))
            f_peptide = str(cmd.get_fastastr("p_"+file[0]))

            f = os.path.join(path, "receptors", file[0]+"_receptor.fa")
            f2 = os.path.join(path, "peptides", file[0]+"_peptide.fa")
                             
            fa_file = open(f, "w")
            fa_file.write(f_protein)
            fa_file.close()
                
            fa_file = open(f2, "w")
            fa_file.write(f_peptide)
            fa_file.close()

            cmd.save(os.path.join(path, "receptors", file[0]+"_receptor.pdb"),"r_"+file[0],-1)
            cmd.save(os.path.join(path, "peptides", file[0]+"_peptide.pdb"),"p_"+file[0],-1)
    
#################################################################################################################################################	

def calculate_bsa():
    cols = ['pdb_id', 'chain 1', 'chain 2', 'protein_sasa', 'peptide_sasa', 'combined_sasa', 'bsa']
    lst2 = []
    
    for filename in os.listdir(os.path.join(path, "solvated_fixed_isolated_pdbs")):

        if filename.endswith(".pdb") or filename.endswith(".cif"):
#             try:
            filenamefull = os.path.join(path, "solvated_fixed_isolated_pdbs", filename)
            cmd.reinitialize()
            cmd.load(filenamefull)

            # gets PDB ID without extension #
            file = os.path.splitext(filename)
            chains = cmd.get_chains(file[0])
            # gets first two chains after removing duplicate entries #
            lst=[]
            for chain in chains:
                count_ca = cmd.count_atoms(file[0]+' and chain ' + chain + ' and name CA')
                lst.append([chain, count_ca])

            lst.sort(reverse=False,key=takeSecond)

            placeholder = ['ZZ', 'YY', 'XX', 'VV', 'WW', 'FF', 'TT', 'PP', 'QQ']

            i = 0
            for chain in chains:
                cmd.alter("chain "+ chain , "chain='"+placeholder[i]+"'" )
                i = i+1
                
            i = 0
            chains = cmd.get_chains(file[0])
            for chain in chains:
                cmd.alter("chain "+ chain , "chain='"+lst[i][0]+"'" )
                i = i+1

            chains = cmd.get_chains(file[0])
            # create objects for protein, peptide, and protein+peptide pair #
            cmd.create("protein", (file[0] + ' and chain '+chains[0]))
            cmd.create("peptide", (file[0] + ' and chain '+chains[1]))
            cmd.create("combined", (file[0] + ' and chain '+chains[0]+'+'+chains[1]))

            # add missing hydrogens #
            cmd.h_add()

            # ignores solvent for calculations #
            cmd.flag("ignore", "none")
            cmd.flag("ignore", "solvent")

            # sets SASA calculation parameters #
            cmd.set("dot_solvent", 1)
            cmd.set("dot_density", 3)

            print("Calculating areas for file ", file[0])

            # calculates areas #
            protein_area=cmd.get_area("protein")
            peptide_area=cmd.get_area("peptide")
            combined_area=cmd.get_area("combined")

            # add to list #
            lst2.append([file[0], chains[0], chains[1], protein_area, peptide_area, combined_area, ((protein_area + peptide_area) - combined_area)])

#             except:
#                 print("failed at PDB ID", filename) 
    
    df1 = pd.DataFrame(lst2, columns=cols)
    dataset_name = os.path.join(path,"areas.csv")
    df1.to_csv(dataset_name)
                
########################################################################################################################################################

if os.path.isdir(args.input):
    path = args.input
    
elif os.path.isfile(args.input):
    path = os.path.join(os.getcwd(), dirName)
    get_pdb_files()
    
if args.clean:
    pdb_separator()
    process_files()

if args.split:
    split_files()

if args.bsa:
    calculate_bsa()


