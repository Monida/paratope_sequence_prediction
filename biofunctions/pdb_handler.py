class PDBHandler:
    '''
    This class allows us to expand on Bio.PDB from Biopython. 
    '''
    def __init__(self,contacts_dict:dict) -> None:
        self.contacts_dict = contacts_dict

    def get_full_seq_from_contacts(self,contact_seq:list,n=0) -> str:
        '''
        Given a sequence of residues (from the antibody or the antigen)
        that are in contact with their counterpart and the PDB ID, get 
        the full residue sequence. 
        N is an additional parameter that inidicates how many residues
        before and after the sequence to add.  
        '''
        full_seq = None
        return full_seq