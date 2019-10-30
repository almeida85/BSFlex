#!/usr/bin/python
#
#
#       Copyright (C) 2010 Yasser Almeida Hernandez <yasser.almeida@gmail.com>
#
#		This script perform the structural superposition of high identity
#		proteins. The reference protein is a protein with a ligand bounded and
#		the rest are ligand free proteins.
#		It uses the Bio.PDB modules in the BioPython package.
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

'''
Superpose Protein Cluster
'''


from Bio.PDB import *
from Bio.SVDSuperimposer import *

# The package p3d is avaible in www.p3d.fufezan.net
# Fufezan, C and Specht, M. (2009): p3d - Python module for structural bioinformatics, BMC Bioinformatics, 10:258
from p3d.vector import Vector

from numpy import array
from math import degrees
import sys, os, string, copy, time

log_file = open('superpose_clusters.log', 'a+')

residues_size_dict = {
    'ALA' : 5,
    'ARG' : 11,
    'ASN' : 8,
    'ASP' : 8,
    'CYS' : 6,
    'GLN' : 9,
    'GLU' : 9,
    'GLY' : 4,
    'HIS' : 10,
    'ILE' : 8,
    'LEU' : 8,
    'LYS' : 9,
    'MET' : 8,
    'PHE' : 11,
    'PRO' : 7,
    'SER' : 6,
    'THR' : 7,
    'TRP' : 14,
    'TYR' : 12,
    'VAL' : 7,
    }

class ProteinCluster():
    '''High Identity Protein Cluster Class'''

    def __init__(self,cluster_name):
        '''ATRIBUTTES
        o ref_struc = Reference structure object (ligand bounded).
        o lig_free_struc_list = List of ligand-free structures.
        '''
        self.lig_free_struc_list = []
        self.ref_struc, self.lig_free_struc_list = self.read_cluster_file(cluster_name)

    ## Private methods
    def _read_structure(self, pdb_code):
        PDB_DIR = '/home/almeida/Travail/PDB/all'

        gunzip_pdb_file = 'gunzip %s/pdb%s.ent.gz' % (PDB_DIR, pdb_code)
        pdb_file = '%s/pdb%s.ent' % (PDB_DIR, pdb_code)

        gzip_pdb_file = 'gzip %s' % pdb_file

        try:
            gzip_pdb_file = 'gzip %s' % pdb_file
        except:
            print pdb_file,': No such file or directory\n'
            log_file.write(pdb_file + ': No such file or directory\n')
            sys.exit()

        ## Decompressing pdb file...
        os.system(gunzip_pdb_file)

        structure = PDBParser()
        PDB = structure.get_structure(string.upper(pdb_code), pdb_file)

        ## Resolution filter
        resolution = PDB.header['resolution']
        if resolution > 2.5:
            #print '\nThe resolution of',string.upper(pdb_code),' exceed the cutoff.'
            log_file.write('The resolution of' + string.upper(pdb_code)+ ' exceed the cutoff\n')
            os.system(gzip_pdb_file)
            return None
        else:
            ## Recompressing pdb file...
            os.system(gzip_pdb_file)
            return PDB


    def _sort_contact_res(self, res_list):
        '''
        Sort a list of residues object according to its sequence number
        '''
        res_index_list = []
        res_list_ord = []

        for res in res_list:
            res_index_list.append(res.id[1])

        resobj_index_list = map(None, res_list, res_index_list)
        resobj_index_list.sort(lambda a,b: cmp(a[1], b[1]))

        for i in resobj_index_list:
            res_list_ord.append(i[0])

        return res_list_ord

    def _get_coordinates(self, atoms):
        '''
        Get the atoms coordinates
        '''
        atoms_coord_list = []
        for atm in atoms:
            atoms_coord_list.append(atm.get_coord())

        coords_array = array(atoms_coord_list, 'f')
        return coords_array

    def _get_backbone_coordinates(self, res):
        '''
        Get the backbone atoms coordinates
        '''

        res_backbone = [res['C'], res['CA'], res['N'], res['CB']]

        return self._get_coordinates(res_backbone)

    def _get_ref_group_atoms_coords(self, res):
        '''
        Get the reference residue GROUP atoms coordinates.
        '''

        if res.resname == 'ARG':
            group_atoms = [res['CD'],res['NE'],res['CZ'],res['NH1'],res['NH2']]
        elif res.resname == 'ASP':
            group_atoms = [res['CG'],res['OD1'],res['OD2']]
        elif res.resname == 'GLU':
            group_atoms = [res['CD'],res['OE1'],res['OE2']]
        elif res.resname == 'TYR':
            group_atoms = [res['CG'],res['CD1'],res['CD2'],res['CE1'],res['CE2'],res['CZ'],res['OH']]

        return self._get_coordinates(group_atoms)

    def _get_equiv_group_atoms_coords(self, res):
        '''
        Get the equivalent residue group atoms coordinates in right and
        inverse order.
        '''
        group_equiv_inv_atoms = []

        ## Get the atoms coordinates in right order
        right_atoms_coords = self._get_ref_group_atoms_coords(res)

        ## Get the atoms coordinates in inverse order
        if res.resname == 'ARG':
            group_equiv_inv_atoms = [res['CD'],res['NE'],res['CZ'],res['NH2'],res['NH1']]
        elif res.resname == 'ASP':
            group_equiv_inv_atoms = [res['CG'],res['OD2'],res['OD1']]
        elif res.resname == 'GLU':
            group_equiv_inv_atoms = [res['CD'],res['OE2'],res['OE1']]
        elif res.resname == 'TYR':
            group_equiv_inv_atoms = [res['CG'],res['CD2'],res['CD1'],res['CE2'],res['CE1'],res['CZ'],res['OH']]

        return right_atoms_coords, self._get_coordinates(group_equiv_inv_atoms)

    def _minimum(self, a, b):
        if a < b:
            return a
        else:
            return b

    def _superpose_guanidine(self, query_struc, ref_res, equiv_res, right_rmsd, inverse_rmsd):
        '''
        Superpose the reference and the equivalent contacting residue based in
        the guanidine group, in right and inverse order (NH1, NH2). Each super-
        position correspond with a RMSD value of the guanidine from the first
        superposition based in the binding sites residues.
        '''
        ref_right_order = [ref_res['CD'], ref_res['NE'], ref_res['CZ'], ref_res['NH1'], ref_res['NH2']]

        equiv_right_order   = [equiv_res['CD'], equiv_res['NE'], equiv_res['CZ'], equiv_res['NH1'], equiv_res['NH2']]
        equiv_inverse_order = [equiv_res['CD'], equiv_res['NE'], equiv_res['CZ'], equiv_res['NH2'], equiv_res['NH1']]

        super_imposer = Superimposer()

        ## Superposition of the equivalent guanidine right order and residue RMSD
        super_imposer.set_atoms(ref_right_order, equiv_right_order)
        super_imposer.apply(query_struc.get_atoms())
        right_guanidine_rmsd = super_imposer.rms

        ## Superposition of the equivalent guanidine inverse order and residue RMSD
        super_imposer.set_atoms(ref_right_order, equiv_inverse_order)
        super_imposer.apply(query_struc.get_atoms())
        inverse_guanidine_rmsd = super_imposer.rms

        if right_guanidine_rmsd < inverse_guanidine_rmsd:
            print 'Guanidine RMSD \t\t\t\t\t\b\b = %.2f A' % right_rmsd
            return right_rmsd
        elif right_guanidine_rmsd > inverse_guanidine_rmsd:
            print 'Guanidine RMSD \t\t\t\t\t\b\b = %.2f A' % inverse_rmsd
            return inverse_rmsd

    def _save_structure(self, name, structure):
        io = PDBIO()
        io.set_structure(structure)
        pdb_out_filename = "%s_ALIGNED.pdb" % name
        io.save(pdb_out_filename)
        return


    ## Public methods
    def read_cluster_file(self, cluster_name):
        '''
        Read the cluster FASTA file. This read the reference structure (holo)
        and the first three ligand-free equivalent structures (apo). Only
        load the structures with resolution <= 2.5
        '''

        cluster_file = open(cluster_name, "r")
        lig_free_struc_list = []
        lig_free_struc_id_list = []

        structure = PDBParser()

        print 'Reading cluster file...'

        for line in cluster_file:
            if line.startswith('>'):
                line = line.strip()
                line = line[1:]
                line_list = line.split("|")

                if line_list[3] == 'CLUSTER_HEAD':
                    pdb_code = string.lower(line_list[0])
                    ref_struc = self._read_structure(pdb_code)

                    ## Resolution filter
                    if ref_struc == None:
                        print 'LOW RESOLUTION!!: The reference structure (',string.upper(pdb_code),') has a low resolution. Cluster ignored\nAborting program\n'
                        log_file.write('LOW RESOLUTION!!: The reference structure ('+string.upper(pdb_code)+') has a low resolution. Cluster ignored\nAborting program\n')
                        sys.exit()

                elif line_list[3] != 'CLUSTER_HEAD':
                    pdb_code = string.lower(line_list[0])
                    lig_free_struc = self._read_structure(pdb_code)

                    ## Resolution filter
                    if lig_free_struc == None:
                        continue

                    elif lig_free_struc.id in lig_free_struc_id_list:
                        continue
                    else:
                        lig_free_struc_id_list.append(lig_free_struc.id)
                        lig_free_struc_list.append(lig_free_struc)
                        if len(lig_free_struc_list) == 3:
                            return ref_struc, lig_free_struc_list
            else:
                continue

        return ref_struc, lig_free_struc_list

    def read_contacts_tables(self, aa):
        '''
        Read the contacts table and save it in a dictionary
        '''

        try:
            cloudy_contacts_table_name = "cloudy_contacts_table_%s.out" % aa
            cloudy_contacts_table = open(cloudy_contacts_table_name, "r")
            contacts_table = open("total_contacts_sidechain-primary.out", "r")
        except:
            print '\nERROR: Some contacts table missing...!!!\n'
            sys.exit()

        cloudy_contacts_dict = {}
        contacts_dict = {}
        id1 = 1
        id2 = 1

        print "\nReading the cloudy contacts table..."
        for line in cloudy_contacts_table:
            line = line.strip()
            line_list = line.split()
            cloudy_contacts_dict[id1] = ''
            cloudy_contacts_dict[id1] = line_list
            id1 += 1

        print "Reading the all contacts table..."
        for line in contacts_table:
            line = line.strip()
            line_list = line.split()
            contacts_dict[id2] = ''
            contacts_dict[id2] = line_list
            id2 += 1

        print 'Cloudy contacts:', len(cloudy_contacts_dict.keys())
        print 'Total contacts :', len(contacts_dict.keys())

        cloudy_contacts_table.close()
        contacts_table.close()

        return cloudy_contacts_dict, contacts_dict

    def check_residue_integrity(self, res):
        '''
        Check for incompletes residues
        '''
        if len(res.get_list()) != residues_size_dict[res.resname]:
            return 1
        else:
            return 0

    def remove_hydrogens(self, res):
        '''
        Remove the hydrogens atoms
        '''

        for atom in res.get_list():
            if atom.id[0] == 'H':
                res.detach_child(atom.id)

        return res

    def get_contact_data(self, cloudy_contacts_dict, all_contacts_dict):
        '''
        This method get the ligand-residue atoms pair objects
        and a list of the contacting residues object, all in the
        Biopython format.
        '''

        contacting_res_list = []
        contacting_res_identifiers_list = []
        binding_site_contacts_list = []

        for cloudy_contact in cloudy_contacts_dict.values():
            pdb_id = cloudy_contact[0].split(".")[0].upper()

            ## Get the "cloudy" contact in the reference structure..
            if pdb_id == self.ref_struc.id:

                ## CONTACTS DATA
                ### Ligand fields...
                ligand_atom_number = cloudy_contact[1]
                ligand_atom_name = cloudy_contact[2]
                ligand_name = cloudy_contact[3]
                ligand_number = int(cloudy_contact[4])
                ligand_chain = cloudy_contact[5]
                ligand_atype = cloudy_contact[6]

                ### Protein fields...
                residue_atom_number = cloudy_contact[7]
                residue_atom_name = cloudy_contact[8]
                residue_name = cloudy_contact[9]
                residue_number = int(cloudy_contact[10])
                residue_chain = cloudy_contact[11]

                ### Contact distance between ligand and residue atoms...
                contact_distance = cloudy_contact[12]

                ## Setting the contacting ligand-residue atoms pair data in the Biopython format...
                ### Residue Data...
                residue_contacting_atom = self.ref_struc[0][residue_chain][residue_number][residue_atom_name]
                residue_contacting = self.ref_struc[0][residue_chain][residue_number]
                residue_contacting_chain = self.ref_struc[0][residue_chain]

                ### Ligand Data...
                cligand_chain = self.ref_struc[0][ligand_chain]
                cligand_name = "H_%s" % ligand_name
                cligand_tuple = (cligand_name,ligand_number," ")

                if self.ref_struc[0][ligand_chain].has_id(cligand_tuple) == 1:
                    cligand_id = self.ref_struc[0][ligand_chain][cligand_tuple]
                    ligand_contacting_atom = self.ref_struc[0][ligand_chain][cligand_tuple][ligand_atom_name]
                else:
                    cligand_id = self.ref_struc[0][ligand_chain][ligand_number]
                    ligand_contacting_atom = self.ref_struc[0][ligand_chain][ligand_number][ligand_atom_name]
            else:
                continue

        ## Get all the contacts in the reference structure..
        for bs_contact in all_contacts_dict.values():
            '''Fields in the contacts table...'''
            pdb_id = bs_contact[0].split(".")[0].upper()

            if pdb_id == self.ref_struc.id:

                ## CONTACTS DATA
                ### Protein fields...
                bs_residue_atom_number = bs_contact[7]
                bs_residue_atom_name = bs_contact[8]
                bs_residue_name = bs_contact[9]
                bs_residue_number = int(bs_contact[10])
                bs_residue_chain = bs_contact[11]

                ### Contact distance between ligand and residue atoms...
                bs_contact_distance = bs_contact[12]

                ## Setting the binding sites contacting residues data in the Biopython format...
                ### Residue Data...
                bs_residue_contacting = self.ref_struc[0][bs_residue_chain][bs_residue_number]

                binding_site_contacts_list.append(bs_residue_contacting)
            else:
                continue

        return residue_contacting_atom, residue_contacting, residue_contacting_chain, cligand_id, ligand_contacting_atom, self._sort_contact_res(binding_site_contacts_list)

    def search_equiv_cont_res(self, ResAtmCont, ResCont, ChainCont, LigID, LigAtmCont, ResContIDs_list):
        '''
        This method get contacting residues list for select this residues in
        the ligand-free structures. First, it check for the presence of this
        residues and then check the number of atoms.
        '''

        FreeProteinsHypContRes_dict = {}  # Conserved contacting residues in the ligand-free proteins

        ## Removing residues wich chain is diferent to contacting chain
        for res in ResContIDs_list:
            if res.get_parent() != ChainCont:
                ResContIDs_list.remove(res)
            else:
                continue

        ## Search for contacting residues in the ligand-free structures
        for struc in self.lig_free_struc_list:
            print '==>',struc.id,'<=='
            log_file.write('\n==> ' + struc.id + ' <==\n')

            for model in struc:
                for chain in model:
                    ## Check if the position of the principal contacting residue exist...
                    if not chain.has_id(ResCont.id[1]):
                        print 'Chain',chain.id,': ERROR!! THE POSITION',ResCont.id[1],'IN',struc.id,'CHAIN',chain.id,'DOESN\'T EXIST. CHAIN IGNORED.'
                        log_file.write('Chain '+chain.id+' : ERROR!! THE POSITION '+str(ResCont.id[1])+' IN '+struc.id+' CHAIN '+chain.id+' DOESN\'T EXIST. CHAIN IGNORED.\n')
                        continue

                    ## Check if the principal contacting residue name is conserved...
                    if chain[ResCont.id[1]].resname != ResCont.resname:
                        print 'Chain',chain.id,': ERROR!! THE EQUIVALENT CONTACTING RESIDUE NAME IS DIFFERENT. CHAIN IGNORED.',\
                        ResCont.resname, ResCont.id[1],'<-->',chain[ResCont.id[1]].resname, ResCont.id[1]
                        log_file.write('Chain '+chain.id+' : ERROR!! THE EQUIVALENT CONTACTING RESIDUE NAME IS DIFFERENT; CHAIN IGNORED. ')
                        log_file.write('('+ResCont.resname+' '+' '+str(ResCont.id[1])+' <--> '+chain[ResCont.id[1]].resname+' '+' '+str(ResCont.id[1]) + ')\n')
                        continue

                    ## Check if the principal contacting residue is complete...
                    if self.check_residue_integrity(chain[ResCont.id[1]]) == 1:
                        print 'Chain',chain.id,' : FATAL ERROR!!: THE RESIDUE IS INCOMPLETE'
                        log_file.write('Chain ' + chain.id + ' : FATAL ERROR!!: THE RESIDUE IS INCOMPLETE\n')
                        continue

                    chain_key = '%s_%s' % (struc.id, chain.id)
                    FreeProteinsHypContRes_dict[chain_key] = ""
                    cont_res_eqiv_backbone = []


                    ## Removing hydrogens and ignoring solvent and ligands...
                    residues = []
                    for res in chain:
                        self.remove_hydrogens(res)
                        if res.id[0][0] == "H" or res.id[0][0] == "W":
                            continue
                        else:
                            residues.append(res)

                    ## Looping by all the residues...
                    nresidues = len(residues)

                    for res in residues:
                        for cont_res in ResContIDs_list:
                            if cont_res.id[1] == res.id[1]:
                                if cont_res.resname == res.resname:
                                    if self.check_residue_integrity(res) == 1:
                                        print 'FATAL ERROR!!: THE CONTACTING RESIDUE IS INCOMPLETE\nABORTING PROGRAM.'
                                        log_file.write('FATAL ERROR!!: THE CONTACTING RESIDUE IS INCOMPLETED\nABORTING PROGRAM.\n')
                                        continue
                                    else:
                                        if res.resname == 'GLY':
                                            cont_res_eqiv_backbone.append(res['C'])
                                            cont_res_eqiv_backbone.append(res['CA'])
                                            cont_res_eqiv_backbone.append(res['N'])
                                        else:
                                            cont_res_eqiv_backbone.append(res['C'])
                                            cont_res_eqiv_backbone.append(res['CA'])
                                            cont_res_eqiv_backbone.append(res['N'])
                                            cont_res_eqiv_backbone.append(res['CB'])

                                else:
                                    if self.check_residue_integrity(res) == 1:
                                        print 'FATAL ERROR!!: THE CONTACTING RESIDUE IS INCOMPLETE\nABORTING PROGRAM.'
                                        log_file.write('FATAL ERROR!!: THE CONTACTING RESIDUE IS INCOMPLETED\nABORTING PROGRAM.\n')
                                        continue
                                    else:
                                        print 'Chain',chain.id,': WARNING: The',cont_res.resname, res.id[1],'IS MUTATED:',cont_res.resname,'<-->',res.resname
                                        log_file.write('Chain '+chain.id+': WARNING: The '+cont_res.resname+' '+str(res.id[1])+' IS MUTATED: '+cont_res.resname+' <--> '+res.resname +'\n')

                                        if res.resname == 'GLY':
                                            cont_res_eqiv_backbone.append(res['C'])
                                            cont_res_eqiv_backbone.append(res['CA'])
                                            cont_res_eqiv_backbone.append(res['N'])
                                        else:
                                            cont_res_eqiv_backbone.append(res['C'])
                                            cont_res_eqiv_backbone.append(res['CA'])
                                            cont_res_eqiv_backbone.append(res['N'])
                                            cont_res_eqiv_backbone.append(res['CB'])
                            else:
                                continue

                    FreeProteinsHypContRes_dict[chain_key] = cont_res_eqiv_backbone
            print

        return FreeProteinsHypContRes_dict

    def superposeAndsave(self, ResAtmCont, ResCont, ResContIDs_list, FreeProteinsHypContResSelect_dict, ChainCont, LigAtmCont):
        '''
        This method superpose the ligand-free proteins on the reference
        complexed protein, based on the selection of the backbone atoms
        of contacting residues in the reference structure. It compute:
        - binding site RMSD.
        - residues chemical group RMSD.
        - contacting residues backbone RMSD.
        - CA-CB(reference)/CA'-CB'(alternate) vectors angle.
        It also superpose the reference and the equivalent contacting residues
        side chains. The experimental contact is saved with each contacting
        residue, so when the last is moved in the superposition, the
        contact is also moved.
        '''

        ## Save the contacting residue with the ligand contacting atom
        #self.get_contacting_ref_resid(self.ref_struc, LigAtmCont, ResCont, ChainCont)

        ## Results file...
        rmsd_file_name = 'clusters_rmsd_%s.out' % ResCont.resname
        rmsd_file = open(rmsd_file_name, 'a+')

        ResCont_backbone = []
        for res in ResContIDs_list:
            if res.resname == 'GLY':
                ResCont_backbone.append(res['C'])
                ResCont_backbone.append(res['CA'])
                ResCont_backbone.append(res['N'])
            else:
                ResCont_backbone.append(res['C'])
                ResCont_backbone.append(res['CA'])
                ResCont_backbone.append(res['N'])
                ResCont_backbone.append(res['CB'])

        log_file.write('\n')

        nchain = 0
        for chain in FreeProteinsHypContResSelect_dict.keys():
            nchain += 1
            if nchain > 3:
                break
            hyp_cont_res_backbone = FreeProteinsHypContResSelect_dict[chain]
            struc_id = chain[:4]
            chain_id = chain[-1]

            ## Superposing structures according to all contacting residues...
            print '==> SUPERPOSING',chain,'ON %s_%s' % (self.ref_struc.id, ChainCont.id)
            log_file.write('==> SUPERPOSING '+chain+' ON '+ self.ref_struc.id+'_'+ChainCont.id+'\n')

            if len(ResCont_backbone) != len(hyp_cont_res_backbone):
                print 'ERROR!!: Fixed and moving atom lists differ in size.\n'
                log_file.write('ERROR!!: Fixed and moving atom lists differ in size. Chain ignored...\n')
                continue

            super_imposer = Superimposer()
            super_imposer.set_atoms(ResCont_backbone, hyp_cont_res_backbone)

            for query_struc in self.lig_free_struc_list:
                if query_struc.id == struc_id:

                    ## Reference vs. Equivalent Binding Sites RMSD
                    super_imposer.apply(query_struc.get_atoms())
                    binding_site_rmsd = super_imposer.rms
                    print 'Binding site contacting residues backbone RMSD = %.2f A' % binding_site_rmsd

                    ## Save aligned structures
                    #self._save_structure(chain, query_struc)

                    ## Compute the RMSD between the contacting residue and the equivalent in the ligand-
                    ## free proteins.
                    ## 'EquivResCont' : Equivalent principal contacting residue in ligand-free proteins
                    ## 'ref_cont_res_coords'   : Reference contacting residue coordinates in the
                    ##                           complexed protein
                    ## 'EquivResCont_coords' : Equivalent principal contacting residue coordinates in
                    ##                           ligand-free proteins
                    ## 'ref_cont_res_backbone, EquivResCont_backbone' : Reference and equivalent
                    ##                                                    backbone residue coordinates
                    EquivResCont = query_struc[0][chain[-1]][ResCont.id[1]]

                    ref_res_group_atoms_coords = self._get_ref_group_atoms_coords(ResCont)
                    equiv_right_group_atoms_coords, equiv_inv_group_atoms_coords = self._get_equiv_group_atoms_coords(EquivResCont)

                    ref_cont_res_backbone = self._get_backbone_coordinates(ResCont)
                    EquivResCont_backbone = self._get_backbone_coordinates(EquivResCont)

                    ## REFERENCE vs. EQUIVALENT RMSD GROUP RMSD...
                    group_rmsd = SVDSuperimposer()

                    ## Right order RMSD
                    group_rmsd.set(ref_res_group_atoms_coords, equiv_right_group_atoms_coords)
                    right_group_rmsd = group_rmsd.get_init_rms()
                    #print "Right order group RMSD\t\t\t\t\b= %.2f A" % right_group_rmsd

                    ## Inverse order RMSD
                    group_rmsd.set(ref_res_group_atoms_coords, equiv_inv_group_atoms_coords)
                    inv_group_rmsd = group_rmsd.get_init_rms()
                    #print "Inverse order group RMSD \t\t\t\b= %.2f A" % inv_group_rmsd

                    if ResCont.resname == 'ARG':
                        ## Superpose reference and equivalent guanidine group and get the correspondent
                        ## right RMSD
                        correct_group_rmsd = self._superpose_guanidine(query_struc, ResCont, EquivResCont, right_group_rmsd, inv_group_rmsd)
                    else:
                        correct_group_rmsd = self._minimum(right_group_rmsd, inv_group_rmsd)

                        print "Minimal group RMSD\t\t\t\t\b= %.2f A" % correct_group_rmsd


                    ## Reference vs. Equivalent Residues Backbones RMSD
                    cont_res_backbone_rmsd = SVDSuperimposer()
                    cont_res_backbone_rmsd.set(ref_cont_res_backbone, EquivResCont_backbone)
                    residues_backbone_rmsd = cont_res_backbone_rmsd.get_init_rms()
                    print "Contacting residues backbone RMSD \t\t\b= %.2f A" % residues_backbone_rmsd

                    ## Angle between the CA-CB vectors in the reference and alternate structures
                    ref_equiv_CACB_vectors_angle = self.compute_CB_CA_vectors_angle(ResCont, EquivResCont)

                    ## Distance between the reference/equivalent protein's residue contacting atom and ligand atom...
                    #exp_distance = ResAtmCont - LigAtmCont
                    #theor_distance = EquivResCont[ResAtmCont.id] - LigAtmCont

                    ## Write the the equivalent contacting residue with the experimental contact
                    #self.get_contacting_equiv_resid(self.ref_struc, query_struc, LigAtmCont, ResCont, EquivResCont, chain[-1], ChainCont)

                    ## Write the RMSD, angles and distances values
                    rmsd_values_string = '%s\t%s\t%s\t%s\t%s\t%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f' % (self.ref_struc.id, ChainCont.id, query_struc.id,\
                    chain_id, ResCont.resname, ResCont.id[1], len(ResContIDs_list), binding_site_rmsd, correct_group_rmsd, residues_backbone_rmsd,\
                    ref_equiv_CACB_vectors_angle)

                    rmsd_file.write(rmsd_values_string + '\n')

        log_file.write('\n\n')
        rmsd_file.close()
        return

    def compute_CB_CA_vectors_angle(self, ResCont, EquivResCont):
        '''
        This method compute the angle between the CA-CB vectors in the reference
        and alternate structure.
        '''

        ref_CA = Vector(ResCont.get_list()[1].coord[0], ResCont.get_list()[1].coord[1], ResCont.get_list()[1].coord[2])
        ref_CB = Vector(ResCont.get_list()[4].coord[0], ResCont.get_list()[4].coord[1], ResCont.get_list()[4].coord[2])

        alt_CA = Vector(EquivResCont.get_list()[1].coord[0], EquivResCont.get_list()[1].coord[1], EquivResCont.get_list()[1].coord[2])
        alt_CB = Vector(EquivResCont.get_list()[4].coord[0], EquivResCont.get_list()[4].coord[1], EquivResCont.get_list()[4].coord[2])

        ref_vector_module = ref_CB.__sub__(ref_CA)
        alt_vector_module = alt_CB.__sub__(alt_CA)

        ref_alt_vector_angle = degrees(ref_vector_module.angleBetween(alt_vector_module))

        return ref_alt_vector_angle

    def get_contacting_ref_resid(self, ref_struc, contact, ref_res, chain_cont):
        '''
        Write the reference contacting residue and the contacting ligand atom
        as a single pdb file.
        '''

        class RefResSelector(Select):
            def accept_residue(self, residue):
                if residue.resname == RefRes.resname and residue.id[1] == RefRes.id[1] and residue.get_parent().id == RefRes.get_parent().id:
                    return 1

        class ExpContact(Select):
            def accept_atom(self,atom):
                if atom.id == contact.id and atom.get_serial_number() == contact.get_serial_number():
                    return 1

        RefStruc = copy.deepcopy(ref_struc)
        RefChainCont = RefStruc[0][chain_cont.id]
        RefRes = RefChainCont[ref_res.id[1]]

        ref_output_name = '%s_%s_reference.pdb' % (RefStruc.id, RefChainCont.id)

        RefResStruc = RefResSelector()

        io=PDBIO()
        io.set_structure(RefStruc)
        io.save(ref_output_name, RefResSelector())
        io.save(ref_output_name, ExpContact())

        return RefRes

    def get_contacting_equiv_resid(self, ref_struc, query_struc, contact, ref_res, EquivResCont, equiv_chain_cont, chain_cont):
        '''
        Write the alternate contacting residue and the contacting ligand atom
        as a single pdb file.
        '''

        class EquivResSelector(Select):
            def accept_residue(self, residue):
                return residue.resname == EquivRes.resname and residue.id[1] == EquivRes.id[1] and residue.get_parent().id == EquivRes.get_parent().id

        class ExpContact(Select):
            def accept_atom(self,atom):
                if atom.id == contact.id and atom.get_serial_number() == contact.get_serial_number():
                    return 1

        EquivStruc = copy.deepcopy(query_struc)
        EquivChainCont = EquivStruc[0][equiv_chain_cont]
        EquivRes = EquivChainCont[ref_res.id[1]]

        alt_output_name = '%s_%s_on_%s_%s.pdb' % (EquivStruc.id, EquivChainCont.id, self.ref_struc.id, chain_cont.id)

        io=PDBIO()
        io.set_structure(EquivStruc)
        io.save(alt_output_name, EquivResSelector())

        io.set_structure(self.ref_struc)
        io.save(alt_output_name, ExpContact())

        return EquivRes


# Functions
def usage():
    '''
    HELP
    '''
    print '\nUSAGE: %s <cluster_fasta_file.fasta> <aa>' % sys.argv[0]
    print '<cluster_fasta_file.fasta>: Cluster file in FASTA format.'
    print '<aa>: Aminoacid to be analyzed (in uppercase). Allowed Only ARG, ASP, GLU and TYR.\n'

def run():
    '''
    RUN SCRIPT
    '''

    if len(sys.argv) != 3:
        usage()
        sys.exit()
    else:
        cluster_name = sys.argv[1]
        allowed_aa = ['ARG', 'ASP', 'GLU', 'TYR']

        if sys.argv[2] not in allowed_aa:
            print 'ERROR!!: Wrong aminoacid.'
            log_file.write('ERROR!!: Wrong aminoacid.\n')
            usage()
            sys.exit()
        else:
            aa = sys.argv[2]

        ## Make the protein cluster object
        Cluster = ProteinCluster(cluster_name)

        ## Store the contacts information in a dictionary
        cloudy_contacts_dict, all_contacts_dict = Cluster.read_contacts_tables(aa)


        ## Get the contact data (as Biopython's objects):
        ## 'ResAtmCont'     : Contacting residue atom
        ## 'ResCont'        : Contacting residue
        ## 'ChainCont'      : Contacting chain
        ## 'LigID'          : Ligand identifier
        ## 'LigAtmCont'     : Contacting ligand atom
        ## 'ResContIDs_list': List with the contacting residues objects
        ResAtmCont, ResCont, ChainCont, LigID, LigAtmCont, ResContIDs_list = Cluster.get_contact_data(cloudy_contacts_dict, all_contacts_dict)

        print '\n======> CLUSTER',Cluster.ref_struc.id,'-',ChainCont.id,'<======'
        print '------- LIGAND DATA --------'
        print 'Ligand name:', LigID.resname,'\nContacting atom:',LigAtmCont.id,'\n----------------------------\n'
        print '------- PROTEIN DATA -------'
        print 'Contacting atom:',ResAtmCont.id,'\nContacting residue:',ResCont.resname,ResCont.id[1],'\nContacting chain:',ChainCont.id
        print '----------------------------\n'

        log_file.write('======> CLUSTER ' + Cluster.ref_struc.id + ' - ' + ChainCont.id + ' <======\n')
        log_file.write('------- LIGAND DATA ---------\n')
        log_file.write('Ligand name: ' + LigID.resname +'\nContacting atom: ' + LigAtmCont.id + '\n-----------------------------\n')
        log_file.write('------- PROTEIN DATA --------\n')
        log_file.write('Contacting atom: ' + ResAtmCont.id + '\nContacting residue: ' + ResCont.resname +' '+' '+ str(ResCont.id[1]) + '\nContacting chain: ' + ChainCont.id+'\n')
        log_file.write('-----------------------------\n')

        ## Dictionary with the equivalent backbone contacting atoms in the ligand-free proteins
        FreeProteinsHypContRes_dict = Cluster.search_equiv_cont_res(ResAtmCont, ResCont, ChainCont, LigID, LigAtmCont, ResContIDs_list)

        ## Superpose the structures, and save the translated coordinates
        Cluster.superposeAndsave(ResAtmCont, ResCont, ResContIDs_list, FreeProteinsHypContRes_dict, ChainCont, LigAtmCont)

        return 0

# MAIN
if __name__ == "__main__":
    run()
    log_file.close()
