#Gets entries for genes of interest from KEGG API, extracts information and generates a tsv file with results

import requests
import gzip
import pickle
import csv

#Creating a Class for each gene and with attributes of interest
class Gene_Classification:
    def __init__(self,gene_name, uniprot_ID, KEGG_ID= '', HRD= False,MSI= False, KEGG_entry=''):
        self.name = gene_name
        self.uniprot_ID= uniprot_ID
        self.HRD = HRD
        self.MSI = MSI
        self.KEGG_ID = KEGG_ID
        self.KEGG_entry= KEGG_entry


print('Extracting Gene Set from gz file...')
#Extracting gene name and uniprot_ID from a Uniprot formatted file
genes_dict = {} #key= gene name, value= uniprot_ID
genes_uniprot_entry_file = gzip.open("gene_set_uniprot_IDs.gz", 'rb')
genes_uniprot_entry_file.readline() #skip header
for line in genes_uniprot_entry_file:
    my_line = line.decode('utf-8') #converting from bytes to string
    entry = str(my_line).split('\t') #creating list out of row
    gene_name = entry[0]
    uniprot_ID = entry[1]
    genes_dict[gene_name]= uniprot_ID


#Create an object for each gene and store it in a list
list_of_objects = []
for gene_name,uniprot_ID in genes_dict.items():
    list_of_objects.append(Gene_Classification(gene_name,uniprot_ID))


#Get KEGG_ID with converting tool from KEGG API
print('Converting Uniprot IDs to KEGG IDs using KEGG API...')
for gene in list_of_objects:
    ID = gene.uniprot_ID
    response = requests.get('http://rest.kegg.jp/conv/genes/uniprot:'+ID)
    try:
        kegg_id = str(response.text).split()[1]
        gene.KEGG_ID = kegg_id
    except: #unknown uniprot IDs
        gene.KEGG_ID = 'not found'
print('Conversion complete.')


print('Getting Entry for each gene from KEGG API...')
#Retrieve entry about each gene and save it as an attribute
#check if defined criteria exists in entry to set attributes

for gene in list_of_objects:
    gene.KEGG_entry= requests.get('http://rest.kegg.jp/get/'+gene.KEGG_ID).text
    if "Mismatch repair" in gene.KEGG_entry:
        gene.MSI = True
    if "Homologous recombination" in gene.KEGG_entry:
        gene.HRD = True
    if "Fanconi anemia" in gene.KEGG_entry:
        gene.HRD = True
        gene.MSI = True


#Saving list of objects in textfile
with open("list_of_gene_objects.txt", "wb") as fp:   #Pickling
    pickle.dump(list_of_objects, fp)


print('Writing results in TSV file...')
with open('KEGG_entry_analysis_results.tsv', 'w') as output_file:
    tsv_writer = csv.writer(output_file, delimiter='\t')
    tsv_writer.writerow(['gene','Uniprot_ID','KEGG_ID', 'MSI', 'HRD'])
    for gene in list_of_objects:
        tsv_writer.writerow([gene.name,gene.uniprot_ID,gene.KEGG_ID,gene.MSI, gene.HRD])