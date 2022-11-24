# paratope_sequence_prediction
Given an antigen predict the aminoacids the paratope needs to have to detect the epitope 

# 1. Database

# 2. Download dataset
[Get the dataset of antibody-antigen complexes with curated affinity data?](http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/about/#faq)

# 3. Understanding the files
http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html

https://github.com/nglviewer/nglview

# 4. Create training dataset

## 4.1 Find epitope and paratope

1. Extract polipeptides
2. Pair-wise comparison of all aminoacids in the antigen and in the interacting AB chain.
3. Sort them by lowest distance first
4. (An epitoge is in general 5 or 6 aminoacids in length)[https://www.pacificimmunology.com/resources/antibody-introduction/what-is-an-epitope/], so we keep the 8 aminoacids with closest proximity to the AB.

## 4.2 Create input-output dataset


# 5. Training the model
https://discuss.huggingface.co/t/trainer-vs-seq2seqtrainer/3145/2
We will use the Trainer class


# Related articles
(Antibody and Antigen Contact Residues Define Epitope and Paratope Size and Structure)
[https://www.jimmunol.org/content/191/3/1428#:~:text=The%20most%20frequent%20size%20of,sizes%20reported%20in%20the%20literature.]