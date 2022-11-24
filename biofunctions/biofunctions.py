"""
This module contains functions to deal specifically with biological data.
"""
import pandas as pd
import torch
from Bio.PDB import Model
from Bio.PDB import PDBParser
import nglview as nv
import datasets


def seqid_is_cdr(ab_res_seqid:int)->tuple:
    """Returns whether an antibody residue sequence ID is in a CDR ,
    as well as the CDR , or an empty string if it is not from
    a CDR.

    Parameters
    ----------
    ab_res_seqid : IMGT sequence identifier (obtained from the get_id() 
    function of the residue object in Bio.PDB module of BioPython.)

    Output
    ------
    is_cdr : bool
        True if the antibody residue sequence ID is in a CDR  and False 
    otherwise.
    cdr : str
        The CDR region the ab_res_seqid belongs to ('CDR1','CDR2','CDR3') or 
    an empty string if it doesn't belong to any CDR region.
    """

    cdr1_start = 27
    cdr1_end = 38
    cdr2_start = 56
    cdr2_end = 65
    cdr3_start = 105
    cdr3_end = 117

    if cdr1_start <= ab_res_seqid and ab_res_seqid <= cdr1_end:
        is_cdr = True
        cdr = 'CDR1'

    elif cdr2_start <= ab_res_seqid and ab_res_seqid <= cdr2_end:
        is_cdr = True
        cdr = 'CDR2'
        
    elif cdr3_start <= ab_res_seqid and ab_res_seqid <= cdr3_end:
        is_cdr = True
        cdr = 'CDR3'

    else:
        is_cdr = False
        cdr = ''

    return is_cdr, cdr


def _add_interaction(ab_label:str, ab_res_name:str, ab_seqid:str, 
        ab_letter:str, ab_letter_seqid:str, cdr:str, ab_list_idx:str,
        ag_label:str, ag_res_name:str, ag_seqid:str, ag_letter:str, 
        ag_letter_seqid:str, ag_list_idx:str,dist:float)->None:
    # Adds the Ab-Ag pair to the globl interactions_dict and the rest of the 
    # passed information to the corresponding global lists.
    # 
    # Paramteres:
    # -----------
    # To understand the paramteres the read get_closes_ress() docstring.
    global compar_count
    ab_labels.append(ab_label)
    ab_ress.append(ab_res_name)
    ab_seqids.append(ab_seqid)
    ab_letters.append(ab_letter)
    ab_letter_seqids.append(ab_letter_seqid)
    cdrs.append(cdr)
    ab_list_idxs.append(ab_list_idx)
    ag_labels.append(ag_label)
    ag_ress.append(ag_res_name)
    ag_seqids.append(ag_seqid)
    ag_letters.append(ag_letter)
    ag_letter_seqids.append(ag_letter_seqid)
    ag_list_idxs.append(ag_list_idx)
    distances.append(dist)
    compar_count += 1

    if compar_count % 100 == 0:
        print(f'{compar_count} comparisons made so far.')
        
    if ab_letter_seqid not in (interactions_dict['interactions']
            [ab_label][cdr]):

        (interactions_dict['interactions'][ab_label]
            [cdr][ab_letter_seqid]) = {}

    if ag_label not in (interactions_dict['interactions'][ab_label][cdr]
            [ab_letter_seqid]):
        (interactions_dict['interactions'][ab_label][cdr]
                [ab_letter_seqid][ag_label]) = []

    (interactions_dict['interactions'][ab_label][cdr]
            [ab_letter_seqid][ag_label]).append(ag_letter_seqid)
    
    if ag_label not in interactions_dict['ag_list_idx'][ab_label][cdr]:
        (interactions_dict['ag_list_idx'][ab_label][cdr]
                [ag_label]) = {'start_list_idx':-1,
                        'end_list_idx':-1}
        
    if (interactions_dict['ag_list_idx'][ab_label][cdr][ag_label]
            ['start_list_idx']) == -1:
        (interactions_dict['ag_list_idx'][ab_label][cdr]
                [ag_label]['start_list_idx']) = ag_list_idx

    else:
        if ag_list_idx < (interactions_dict['ag_list_idx'][ab_label][cdr]
                [ag_label]['start_list_idx']):
            (interactions_dict['ag_list_idx'][ab_label][cdr]
                    [ag_label]['start_list_idx']) = ag_list_idx
    
    if ag_list_idx > (interactions_dict['ag_list_idx'][ab_label][cdr]
            [ag_label]['end_list_idx']):
        (interactions_dict['ag_list_idx'][ab_label][cdr][ag_label]
                ['end_list_idx']) = ag_list_idx


def get_closest_ress(model:Model, Ab_labels:list, Ag_labels:list, 
        max_dist:int, amino_acids_dict:dict)->tuple:
    """Given a Bio.PDB.Model, compares all the CA atoms of every Ab and
    Ag chain ressidues and keeps the closest Ab-Ag residue pairs (e.g., the
    pairs whose distance is less or equal to max_dist).

    Returns 2 dictionaries: closest_ress_dict and interactions_dict.

    Paramters
    ---------
    Ab_labels : list of Ab chain labels to compare
    Ag_labels : list of Ag chain labels to compare
    max_dist : Max distance to keep an Ab-Ag residure pair.

    Outuput
    -------
    closest_ress_dict: A dictionary where the values are lists, and where each
    element of each list corresponds to a single Ab-Ag residue pair. The 
    hierarchy is as follows:
        key:value pairs
        * 'ab_label':ab_labels (list of Ab chain labels)
        * 'ab_res':ab_ress (list of Ab residue names)
        * 'cdr':cdrs (list of the CDR region the Ab
        residue belongs to.)
        * 'list_idx':list_idxs (lst of the list index of the Ab residue
        within the Ab sequence string)
        * 'ab_seqid':ab_seqids (list of Ab residue sequence ids)
        * 'ag_label':ag_labels (list of Ag chain labels)
        * 'ag_res':ag_ress (list of Ag residue names)
        * 'ag_seqid':ag_seqids (of Ag residue sequence ids)
        * 'distances':distances (list Ab-Ag residue distances)

    interactions_dict: A dictionary that contains more information about each 
    Ab-Ag residue pair that is believed to be in the paratope-epitope 
    interaction. The hierarchy is as follows:
        * full_ab_seq: a dictionary key
            * Ab chain label: the letter corresponding to the Ab chain
                * The Ab chain sequence until the residue with sequence id 
                (seqid) = 117 in the IMGT numbering system. This is chosen as 
                the limit becuase it is the last residue that is considered 
                to be in a CDR, therefore the first 117 seqid are the only 
                residues of interest.
        * full_ab_seq_map: a dictionary key
            * Ab chain label: the letter corresponding to the Ab chain
                * The same Ab chain sequence as in **full_ab_seq**, but every 
                residue is mapped to its corresponding seqid.
        * full_ag_seq: a dictionary key
            * Ag chain label: the letter corresponding to the Ag chain
                * The full Ag sequence of what is considered to be the 
                epitope. 
        * full_ag_seq_map: a dictionary key
            * Ag chain label: the letter corresponding to the Ag chain
                * The same Ag sequence as in **full_ag_seq**, but every 
                residue is mapped to its corresponding seqid.
        * ab_list_idx: a dictionary key
            * Ab chain label: the letter corresponding to the Ab chain
                * CDR: The corresponding CDR ('CDR1', 'CDR2' or 'CDR3')
                    * It is a 2-element tupple that maps to the starting and 
                    ending seqids of the CDR to the position within the 
                    full_ab_seq when read as a python list. The 1st element 
                    maps to the starting seqid and the 2nd element maps to 
                    the ending seqid.
        * ag_list_idx: a dictionary key
            * Ab chain label: the letter corresponding to the Ab chain  
                * CDR: The corresponding CDR ('CDR1', 'CDR2' or 'CDR3')
                    * Ag chain label: the letter corresponding to the Ag 
                    chain that is in contact with the corresponding CDR.
                        * start_list_idx: Maps to the starting seqid of the 
                        Ag sequence.
                        * end_list_idx: Mapst to the ending seqid of the
                        Ag sequence.
        * interactions: a dictionary that 
            * Ab chain label: the letter corresponding to the Ab chain  
                * CDR: The corresponding CDR ('CDR1', 'CDR2' or 'CDR3')
                    * Ab letter-seqid: the letter corresponding to the Ab 
                    residue contatenated to its seqid (e.g., 'Y-27')
                        * Ag chain label:  the letter corresponding to the Ag 
                        chain that is in contact with the corresponding CDR.
                            * Ag letter-seqid: the letter corresponding to 
                            the Ab residue contatenated to its seqid 
                            (e.g., 'Y-27').
    """
    global ab_labels, ab_ress, ab_seqids, ab_letters, ab_letter_seqids
    global cdrs, ab_list_idxs, ag_labels, ag_ress, ag_seqids, ag_letters
    global ag_letter_seqids, ag_list_idxs, distances, compar_count
    global interactions_dict

    ab_labels = []
    ab_ress = []
    ab_seqids = []
    ab_letters = []
    ab_letter_seqids = []
    cdrs = []
    ab_list_idxs = []
    ag_labels = []
    ag_ress = []
    ag_seqids = []
    ag_letters = []
    ag_letter_seqids = []
    ag_list_idxs = []
    distances = []
    compar_count = 0
    interactions_dict = {'full_ab_seq':{},
                    'full_ab_seq_map':{},
                    'full_ag_seq':{},
                    'full_ag_seq_map':{},
                    'ab_list_idx':{},
                    'ag_list_idx':{},
                    'interactions':{}}

    for ab_label in Ab_labels:
        full_ab_seq = ''
        full_ab_seq_map = []
        Ab_ch = model[ab_label]
        ab_list_idx = 0
        cdr = ''
        cdr1_start_idx = None
        cdr1_end_idx = None
        cdr2_start_idx = None
        cdr2_end_idx = None
        cdr3_start_idx = None
        cdr3_end_idx = None

        if ab_label not in interactions_dict['interactions']:
            interactions_dict['interactions'][ab_label] = {'CDR1':{},
                                                    'CDR2':{},
                                                    'CDR3':{}}

        if ab_label not in interactions_dict['ag_list_idx']:
            interactions_dict['ag_list_idx'][ab_label] = {'CDR1':{},
                                                    'CDR2':{},
                                                    'CDR3':{}} 

        for ab_res in Ab_ch.get_residues():
            ab_res_name = ab_res.get_resname()
            ab_res_het_tag, ab_seqid, _ = ab_res.get_id()

            if ab_res_het_tag == ' ' or ab_res_het_tag == '':
                ab_letter = amino_acids_dict[ab_res_name]
                ab_letter_seqid = f'{ab_letter}-{ab_seqid}'
                full_ab_seq = full_ab_seq + ab_letter
                full_ab_seq_map .append(ab_letter_seqid)

                ab_atom = ab_res['CA']

                for ag_label in Ag_labels:
                    Ag_ch = model[ag_label]
                    full_ag_seq = ''
                    full_ag_seq_map = []
                    ag_list_idx = 0
                    
                    for ag_res in Ag_ch.get_residues():
                        ag_res_name = ag_res.get_resname()
                        ag_res_het_tag, ag_seqid, _ = ag_res.get_id()
                        
                        if ag_res_het_tag == ' ' or ag_res_het_tag == '':
                            ag_letter = amino_acids_dict[ag_res_name]
                            ag_letter_seqid = f'{ag_letter}-{ag_seqid}'
                            full_ag_seq = full_ag_seq + ag_letter
                            full_ag_seq_map .append(ag_letter_seqid)

                            ag_atom = ag_res['CA']

                            dist = ab_atom - ag_atom

                            is_cdr, cdr = seqid_is_cdr(ab_seqid)

                            if cdr == 'CDR1':
                                if ab_seqid == 27 and not cdr1_start_idx:
                                    cdr1_start_idx = ab_list_idx
                                if ab_seqid == 38 and not cdr1_end_idx:
                                    cdr1_end_idx = ab_list_idx

                            if cdr == 'CDR2':
                                if ab_seqid == 56 and not cdr2_start_idx:
                                    cdr2_start_idx = ab_list_idx
                                if ab_seqid == 65 and not cdr2_end_idx:
                                    cdr2_end_idx = ab_list_idx

                            if cdr == 'CDR3':
                                if ab_seqid == 105 and not cdr3_start_idx:
                                    cdr3_start_idx = ab_list_idx
                                if ab_seqid == 117 and not cdr3_end_idx:
                                    cdr3_end_idx = ab_list_idx

                            if is_cdr and (ab_label in Ab_labels) and \
                                (dist<=max_dist):

                                _add_interaction(ab_label, ab_res_name,
                                        ab_seqid, ab_letter,
                                        ab_letter_seqid, cdr, ab_list_idx,
                                        ag_label, ag_res_name, ag_seqid,
                                        ag_letter, ag_letter_seqid,
                                        ag_list_idx, dist)

                        ag_list_idx += 1
                        
                    (interactions_dict['full_ag_seq']
                            [ag_label]) = full_ag_seq
                    (interactions_dict['full_ag_seq_map']
                            [ag_label]) = full_ag_seq_map

            ab_list_idx += 1      
              
        if ab_label not in interactions_dict['ab_list_idx']:
            interactions_dict['ab_list_idx'][ab_label] = {}

        interactions_dict['ab_list_idx'][ab_label]['CDR1'] = (cdr1_start_idx, 
                cdr1_end_idx)
        interactions_dict['ab_list_idx'][ab_label]['CDR2'] = (cdr2_start_idx, 
                cdr2_end_idx)
        interactions_dict['ab_list_idx'][ab_label]['CDR3'] = (cdr3_start_idx, 
                cdr3_end_idx)

        interactions_dict['full_ab_seq'][ab_label] = full_ab_seq
        interactions_dict['full_ab_seq_map'][ab_label] = full_ab_seq_map

    print(f'{compar_count} total comparisons')
    print('Finished processing pdb.')

    closest_ress_dict = {'ab_label':ab_labels,
            'ab_res':ab_ress,
            'ab_seqid':ab_seqids,
            'ab_letter':ab_letters,
            'ab_letter_seqid':ab_letter_seqids,
            'cdr':cdrs,
            'ab_list_idx':ab_list_idxs,
            'ag_label':ag_labels,
            'ag_res':ag_ress,
            'ag_seqid':ag_seqids,
            'ag_letter':ag_letters,
            'ag_letter_seqid':ag_letter_seqids,
            'ag_list_idx':ag_list_idxs,
            'distance':distances,
            }

    return closest_ress_dict, interactions_dict


def create_in_out_str(interactions_dict:dict, n:int)->list:
    """For a given interactions_dict, it returns for each interaction,
    a dictionary with the following two keys:
        * input_text:
        * output_text:
    """
    in_out_dict = {}
    in_out_list = []

    for ab_label in interactions_dict['interactions']:

        for cdr in interactions_dict['interactions'][ab_label]:
            cdr_start = interactions_dict['ab_list_idx'][ab_label][cdr][0]
            cdr_end = interactions_dict['ab_list_idx'][ab_label][cdr][1]

            if cdr_start and cdr_end:
                if cdr_start - n < 0:
                    cdr_start = 0
                else:
                    cdr_start = cdr_start - n

                if cdr_end + n > len(interactions_dict['full_ab_seq_map']
                                     [ab_label]):

                    cdr_end = len(interactions_dict['full_ab_seq_map']
                                  [ab_label])

                else:
                    cdr_end = cdr_end + n + 1

                ab_rows = (interactions_dict['full_ab_seq_map'][ab_label]
                           [cdr_start:cdr_end])

                ab_in = (interactions_dict['full_ab_seq'][ab_label]
                        [cdr_start:cdr_end])

                for ag_label in (interactions_dict['ag_list_idx'][ab_label]
                                 [cdr]):

                    ag_start = (interactions_dict['ag_list_idx'][ab_label]
                                [cdr][ag_label]['start_list_idx'])

                    ag_end = (interactions_dict['ag_list_idx'][ab_label]
                              [cdr][ag_label]['end_list_idx'])

                    if ag_start and ag_end:
                        if ag_start - n < 0:
                            ag_start = 0

                        else:
                            ag_start = ag_start - n

                        if ag_end + n > len(interactions_dict['full_ag_seq_map']
                                            [ag_label]):

                            ag_end = len(interactions_dict['full_ag_seq_map']
                                         [ag_label])

                        else:
                            ag_end = ag_end + n + 1

                    ag_cols = (interactions_dict['full_ag_seq_map'][ag_label]
                               [ag_start:ag_end])

                    ag_in = (interactions_dict['full_ag_seq'][ag_label]
                             [ag_start:ag_end])

                    interactions_df = pd.DataFrame(columns=ag_cols\
                                                           +['sentence_end'], 
                                                   index=ab_rows)

                    for ab_letter_seqid, _ in (interactions_dict['interactions']
                                               [ab_label][cdr]).items():

                        if ag_label in (interactions_dict['interactions']
                                        [ab_label][cdr][ab_letter_seqid]):

                            for ag_letter_seqid in (interactions_dict['interactions']
                                                    [ab_label][cdr][ab_letter_seqid]
                                                    [ag_label]):

                                interactions_df.at[ab_letter_seqid,ag_letter_seqid] = '|'

                    interactions_df.fillna('.',inplace=True)
                    interactions_df['sentence_end'] = '_'
                    interactions_arr = interactions_df.to_numpy()
                    size = interactions_arr.size
                    interactions_arr = interactions_arr.reshape(1,size)
                    interactions_str = ''.join(list(interactions_arr[0]))
                    in_out_dict = {'input_text':f'<{cdr}>{ab_in}<ag>{ag_in}',
                                    'output_text': f'<out>{interactions_str}'}
                    in_out_list.append(in_out_dict)

    return in_out_list       


class DataParser:
    def __init__(self, data: datasets.Dataset, vocab:list)->None:
        self.data = data
        self.vocab = vocab
        self.vocab_map = {item:idx for (idx,item) in enumerate(vocab)}
        self.reverse_vocab_map = {idx:item for (idx,item) in enumerate(vocab)}

    def encode(self, s:str)->list:
        """Transform a string of aminoacid letter codes to a list of vocab encoding
        indices.
        """
        l = list(s)
        indices = [self.vocab_map[item] for item in l]
        return indices

    def decode(self, indices:list)->list:
        """Transforms a list of vocab encoding indices into a string of aminoacid
        letter codes.
        """
        labels = [self.reverse_vocab_map[idx] for idx in indices]
        return labels

    def encode_data(self)->None:
        """Encodes all the input_text and output_text of all the data.
        """
        encoded_data = []
        for temp_dict in self.data:
            ag_idx = temp_dict['input_text'].index('<ag>')
            cdr = temp_dict['input_text'][:6]
            input_cdr = temp_dict['input_text'][6:ag_idx]
            input_ag = temp_dict['input_text'][ag_idx+4:]
            input_lab = temp_dict['output_text'][5:]

            temp_dict['input_ids'] = [self.vocab_map[f'{cdr}']] \
                                        + self.encode(input_cdr) \
                                        + [self.vocab_map['<ag>']] \
                                        + self.encode(input_ag)

            temp_dict['labels'] = [self.vocab_map['<out>']] + self.encode(input_lab)
            encoded_data.append(temp_dict)
        self.data = encoded_data
        return encoded_data

    def __getitem__(self, idx:int)->dict:
        return self.data[idx]

    def __len__(self):
        return len(self.data)

class DataCollator:
    def __init__(self, vocab:list)->dict:
        self.vocab = vocab
        self.vocab_map = {item:idx for (idx,item) in enumerate(vocab)}

    def padding(self, in_or_out_list:list)->list:
        """Adds the <bos>, <eos> and <pad> words to a list of input or output
        strings. The padding is added to match the longest input or output.

        Parameters
        ----------
        in_or_out_list : list
            A list where each element is a list of input_indices or label_indices

        Outuput
        -------
        padded_in_or_out_list : list
            A list where each element is a list of padded input or label indices. 
        """
        max_len = max([len(item) for item in in_or_out_list])

        padded_in_or_out_list = [] 

        for item in in_or_out_list:
            pad = max_len - len(item)
            l = [self.vocab_map['<bos>']] + \
                item + \
                [self.vocab_map['<eos>']] + \
                [self.vocab_map['<pad>']]*pad

            padded_in_or_out_list.append(l)

        return padded_in_or_out_list
    
    def attention_masking(self,padded_inputs:list)->list:
        """Adds ones to all entries that represent tokens and zeros to all 
        padding tokens. 
        
        Paramters
        ---------
        padded_inputs : list
            A list where each element is a list of padded input, that where
        the output of self.padding() when a list of input_indices was passed.
        """

        masked_inputs = []
        for item in padded_inputs:
            masked_inputs.append([0 if token == 1 else 1 for token in item])
        return masked_inputs

    def collator_call(self,data_examples_list:list)->dict:
        """For a batch of examples from the data, replace input_ids and labels
        with padded input_ids and labels. It also adds attention mask.

        Parameters
        ----------
        data_examples_list: a list where each element is a dictionary that contains
        the following keys:
            * input_text:
            * output_text:
            * input_ids:
            * labels:

        """
        batch = {}
        input = []
        outputs = []

        # Each element of data_examples_list is a dictionary
        for item in data_examples_list:
            for key in item:
                if key == 'input_ids':
                    input.append(item[key])
                elif key == 'labels':
                    outputs.append(item[key])

        padded_inputs = self.padding(input)
        batch['input_ids'] = torch.tensor(padded_inputs,dtype=torch.long)

        padded_outputs = self.padding(outputs)
        batch['labels'] = torch.tensor(padded_outputs,dtype=torch.long)

        masked_inputs = self.attention_masking(padded_inputs)
        batch['attention_mask'] = torch.tensor(masked_inputs, dtype=torch.long)
        return batch
    
    def __call__(self, data_examples_list:list)->dict:
        return self.collator_call(data_examples_list)

class Tokenizer:
    def __init__(self,vocab:dict):
        self.vocab = vocab
