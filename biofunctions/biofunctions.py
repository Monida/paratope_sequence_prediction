from Bio.PDB import Model

def ab_res_seqid_is_cdr(ab_res_seqid:int)->tuple:
    '''
    Returns whether an antibody residue sequence ID is in a CDR region,
    as well as the CDR region.

    Parameters
    ----------
    ab_res_seqid : IMGT sequence identifier (obtained from the get_id() function of 
    the residue object in Bio.PDB module of BioPython.)

    Output
    ------
    is_cdr : bool
        True if the antibody residue sequence ID is in a CDR region and False otherwise.
    cdr_region : str
        The CDR region the ab_res_seqid belongs to ('CDR1','CDR2','CDR3') or an empty 
    string if it doesn't belong to any CDR region.
    '''
    cdr1_start = 27
    cdr1_end = 38
    cdr2_start = 56
    cdr2_end = 65
    cdr3_start = 105
    cdr3_end = 117

    if cdr1_start <= ab_res_seqid and ab_res_seqid <= cdr1_end:
        is_cdr = True
        cdr_region = 'CDR1'
    elif cdr2_start <= ab_res_seqid and ab_res_seqid <= cdr2_end:
        is_cdr = True
        cdr_region = 'CDR2'
    elif cdr2_start <= ab_res_seqid and ab_res_seqid <= cdr2_end:
        is_cdr = True
        cdr_region = 'CDR3'
    else:
        is_cdr = False
        cdr_region = ''
    return is_cdr, cdr_region

def get_ab_ab_closest_ress(model: Model, Ab_labels:list, Ag_labels:list, max_dist:int)->dict:
    '''
    Given a Bio.PDB.Model, compares all the CA atoms of every Ab and
    Ag ressidue and keeps the closest Ab-Ag residue pairs (e.g., the
    pairs whose distance is less than max_dist).

    Returns a dictionary with the Ab-Ag distances as well as other 
    identifier information.

    Paramters
    ---------
    Ab_labels : list of Ab chain labels to compare
    Ag_labels : list of Ag chain labels to compare
    max_dist : Max distance to keep an Ab-Ag residure pair.

    Outuput
    -------
    comparisons_dict:
        key:value pairs
        * 'ab_label':ab_labels (list of Ab chain labels)
        * 'ab_res':ab_ress (list of Ab residue names)
        * 'ab_seqid':ab_seqids (list of Ab residue sequence ids)
        * 'ag_label':ag_labels (list of Ag chain labels)
        * 'ag_res':ag_ress (list of Ag residue names)
        * 'ag_seqid':ag_seqids (of Ag residue sequence ids)
        * 'distances':distances (list Ab-Ag residue distances)
        Note: The first element of each list corresponds to the first
        Ab-Ag residue comparison, the second element corresponds to the 
        second comparison and so on. 
    '''
    distances = []
    ab_ress = []
    ag_ress = []
    ab_seqids = []
    ag_seqids = []
    ab_labels = []
    ag_labels = []
    cdr_regions = []
    comparison_num = 0
    for ab_label in Ab_labels:
        Ab_ch = model[ab_label]
        for ab_res in Ab_ch.get_residues():
            ab_res_name = ab_res.get_resname()
            ab_res_het_tag, ab_seqid, _ = ab_res.get_id()
            if ab_res_het_tag == ' ' or ab_res_het_tag == '':
                ab_atom = ab_res['CA']
                for ag_label in Ag_labels:
                    Ag_ch = model[ag_label]
                    for ag_res in Ag_ch.get_residues():
                        ag_res_name = ag_res.get_resname()
                        ag_res_het_tag, ag_seqid, _ = ag_res.get_id()
                        if ag_res_het_tag == ' ' or ag_res_het_tag == '':
                            ag_atom = ag_res['CA']
                            dist = ab_atom - ag_atom
                            is_cdr, cdr_region = ab_res_seqid_is_cdr(ab_seqid)
                            if is_cdr and (ab_label in Ab_labels) and (dist<=max_dist):
                                ab_ress.append(ab_res_name)
                                ag_ress.append(ag_res_name)
                                ab_seqids.append(ab_seqid)
                                ag_seqids.append(ag_seqid)
                                ab_labels.append(ab_label)
                                ag_labels.append(ag_label)
                                cdr_regions.append(cdr_region)
                                distances.append(dist) 
                                comparison_num += 1
                                if comparison_num % 100 == 0:
                                    print(f'{comparison_num} comparisons made so far.')


    print(f'{comparison_num} total comparisons')
    print('Finished...')
    comparisons_dict = {
        'ab_label':ab_labels,
        'ab_res':ab_ress,
        'ab_seqid':ab_seqids,
        'ag_label':ag_labels,
        'ag_res':ag_ress,
        'ag_seqid':ag_seqids,
        'distances':distances
        }
    return comparisons_dict