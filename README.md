# paratope_sequence_prediction
Given an antigen predict the aminoacids the paratope needs to have to detect the epitope 

# 1. Database

# 2. Download dataset
[Get the dataset of antibody-antigen complexes with curated affinity data?](http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/about/#faq)

# 3. Understanding the files
http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html

https://github.com/nglviewer/nglview

Added the get_start_end to the Polypeptide class

````
    def get_start_end(self):
        """Return the start and end of the polypeptide.

        Return a tuppple (START,END) where START
        and END are sequence identifiers of the outer residues.
        """
        start = self[0].get_id()[1]
        end = self[-1].get_id()[1]

        return start, end
````

based on 

````
    def __repr__(self):
        """Return string representation of the polypeptide.

        Return <Polypeptide start=START end=END>, where START
        and END are sequence identifiers of the outer residues.
        """
        start = self[0].get_id()[1]
        end = self[-1].get_id()[1]
        return "<Polypeptide start=%s end=%s>" % (start, end)
````
# 4. Find epitope and paratope

## 4.1 Annotate pdb files
1. original pdb to fasta using (Bio.SeqIO)[https://biopython.org/wiki/SeqIO]
````
pdbfile = 'sabdab-data/original_pdb/4liq.pdb'
with open(pdbfile) as handle:
    for record in SeqIO.parse(handle, 'pdb-seqres'):
        print(record.id)
        chain_label = record.id.split(':')[-1]
        with open(f'sabdab-data/fasta/4liq_{chain_label}.fasta', 'w') as output_handle:
            SeqIO.write(record, output_handle, 'fasta')
````

2. Use ANARCI to annotate the chains
http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/anarci/

## 4.2 How to know which chains interact?

1. Extract polipeptides
2. Pair-wise comparison of all aminoacids in the antigen and in the interacting AB chain.
3. Sort them by lowest distance first
4. (An epitoge is in general 5 or 6 aminoacids in length)[https://www.pacificimmunology.com/resources/antibody-introduction/what-is-an-epitope/], so we keep the 8 aminoacids with closest proximity to the AB.

# Related articles
(Antibody and Antigen Contact Residues Define Epitope and Paratope Size and Structure)
[https://www.jimmunol.org/content/191/3/1428#:~:text=The%20most%20frequent%20size%20of,sizes%20reported%20in%20the%20literature.]