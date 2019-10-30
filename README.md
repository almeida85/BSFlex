# Binding sites superposition

This script is part of a study of the induced-fit effect in proteins upon ligand binding. It uses a previously compiled database (BMC Bioinformatics (2007) volume 8, number 310, DOI:10.1186/1471-2105-8-310, and Protein Science (2002), Volume 47, Issue 1, DOI:10.1002/prot.10026), wich stores information about the contacting atoms from proteins residues in the binding site and ligands, from the Protein Data Bank (PDB).

The goal of this study is to characterize the flexibility of the residues in binding sites by asking the following question: how would be the residue/ligand relative geometry if no induced-fit effect occurs?

superpose_cluster.py script perform the structural superposition of high identity proteins identified with BLASTIdentityClust.py (repoSeq respository). It reads the database with the information of the residue-ligand contacts, and extract de coordinates of the original contact ("holo" contact) and the theoretical contact ("apo" contact)

The reference protein is a protein with a ligand bounded and the rest are ligand free proteins.

It uses the Bio.PDB modules in the BioPython package.
