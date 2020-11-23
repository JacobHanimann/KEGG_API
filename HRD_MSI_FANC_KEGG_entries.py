#Gets entries for genes of interest from KEGG API, extracts information and generates a tsv file with results

import requests
import gzip
import pickle
import csv


#Creating a Class for each gene and with attributes of interest
class Gene_Classification:
    def __init__(self,gene_name, uniprot_ID, KEGG_ID= '', HRD= 0,MSI= 0, KEGG_entry='', FANC=0):
        self.name = gene_name
        self.uniprot_ID= uniprot_ID
        self.HRD = HRD
        self.MSI = MSI
        self.FANC = FANC
        self.KEGG_ID = KEGG_ID
        self.KEGG_entry= KEGG_entry


print('Extracting Gene Set from gz file...')
#Extracting gene name and uniprot_ID from a Uniprot formatted file
genes_dict = {} #key= gene name, value= uniprot_ID
genes_uniprot_entry_file = gzip.open("gene_panel_uniprot_IDs.gz", 'rb')
genes_uniprot_entry_file.readline() #skip header
for line in genes_uniprot_entry_file:
    my_line = line.decode('utf-8') #converting from bytes to string
    entry = str(my_line).split('\t') #creating list out of row
    gene_name = entry[0]
    uniprot_ID = entry[1]
    genes_dict[gene_name]= uniprot_ID


#Create an object for each gene and store it in a list
list_of_gene_objects = [Gene_Classification(gene_name, uniprot_ID) for gene_name, uniprot_ID in genes_dict.items()]


#Get KEGG_ID with converting tool from KEGG API
print('Converting Uniprot IDs to KEGG IDs using KEGG API...')
for gene in list_of_gene_objects:
    ID = gene.uniprot_ID
    response = requests.get('http://rest.kegg.jp/conv/genes/uniprot:'+ID)
    try:
        kegg_id = str(response.text).split()[1]
        gene.KEGG_ID = kegg_id
        print(gene.name,gene.KEGG_ID)
    except: #unknown uniprot IDs
        print(gene.name+' not found')
print('Conversion complete.')


print('Getting Entry for each gene from KEGG API...')
#Retrieve entry about each gene and save it as an attribute
#check if defined criteria exists in entry to set attributes

for gene in list_of_gene_objects:
    gene.KEGG_entry= requests.get('http://rest.kegg.jp/get/'+gene.KEGG_ID).text
    if "Mismatch repair" in gene.KEGG_entry:
        gene.MSI = 1
        print('Mismatch gene: '+gene.name)
    if "Homologous recombination" in gene.KEGG_entry:
        gene.HRD = 1
        print('Homologous recombination gene: '+gene.name)
    if "Fanconi anemia" in gene.KEGG_entry:
        gene.FANC=1
        print('Fanconi anemia: '+ gene.name)


#Saving list of objects in textfile
with open("list_of_gene_objects.txt", "wb") as fp:   #Pickling
    pickle.dump(list_of_gene_objects, fp)


print("Complementing data manually...")

for gene in list_of_gene_objects:
    if gene.name =='C11orf30,EMSY':
        gene.HRD =1
    if gene.name == 'ATRX, RAD54L':
        gene.HRD = 1


print('Check genes from studies to complement KEGG data...')
#check if genes of study (Homologous recombination deficiency status-based classification of high-grade serous ovarian carcinoma 2020,https://www.nature.com/articles/s41598-020-59671-3) of HRD are also associated on KEGG
count= 0
undiscovered=0
gene_list_study = ('BRCA1,BRCA2,ATM,ATR,BARD1,BLM,BRIP1,CDK12,CHEK1,CHEK2,FANCA,FANCC,FANCD2,FANCE,FANCF,FANCI,FANCL,FANCM,MRE11,NBN,PALB2,RAD50,RAD51,RAD51B,RAD51C,RAD51D,RAD52,RAD54L,RPA1').split(',')
for gene in list_of_gene_objects:
    if gene.name in gene_list_study:
        count +=1
        if gene.HRD ==0:
            undiscovered +=1
            gene.HRD = 1

print(str(count)+' genes investigated in the study are part of the Cdx gene set and '+str(undiscovered)+' of them were not annotated as associated with Homologous repair in KEGG.')


print('Check genes from studies to complement KEGG data...')
#check if genes of study https://www.nature.com/articles/nrc.2015.21/tables/1 of HRD are also associated on KEGG
count= 0
undiscovered=0
gene_list_study = ('ATM,ATR,BAP1,BRCA1,BRCA2,CDK12,CHEK2,CHK2,FANCA,FANCC,FANCD2,FANCE,FANCF,PALB2,NBS1,NBN,WRN,RAD51C,RAD51D,MRE11A,CHEK1,CHK1,BLM,RAD51B,BRIP1').split(',')
for gene in list_of_gene_objects:
    if gene.name in gene_list_study:
        count +=1
        if gene.HRD ==0:
            undiscovered +=1
            gene.HRD = 1

print(str(count)+' genes investigated in the study are part of the Cdx gene set and '+str(undiscovered)+' of them were not annotated as associated with Homologous repair in KEGG.')


#check if genes of a study (Microsatellite Instability Is Associated With the Presence of Lynch Syndrome Pan-Cancer, 2019, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6553803/) of MSI are also annotated on KEGG
count= 0
undiscovered=0
gene_list_study = ('MLH1,MSH2,MSH6,PMS2,EPCAM').split(',')
for gene in list_of_gene_objects:
    if gene.name in gene_list_study:
        count +=1
        if gene.MSI ==0:
            undiscovered +=1
            gene.MSI = 1

print(str(count)+' genes investigated in the study are part of the Cdx gene set and '+str(undiscovered)+' of them were not annotated as associated with MSI in KEGG.')


print('Overview of pathways:')
HRD = []
MSI = []
FANC = []
for gene in list_of_gene_objects:
    if gene.MSI==1:
        MSI.append(gene.name)
    if gene.HRD==1:
        HRD.append(gene.name)
    if gene.FANC == 1:
        FANC.append(gene.name)
print('HRD: '+str(HRD))
print('MSI: '+str(MSI))
print('FANC: '+str(FANC))


print('Writing results in TSV file...')
with open('KEGG_entry_analysis_results.tsv', 'w') as output_file:
    tsv_writer = csv.writer(output_file, delimiter='\t')
    tsv_writer.writerow(['gene','Uniprot_ID','KEGG_ID', 'MSI', 'HRD','FANC'])
    for gene in list_of_gene_objects:
        tsv_writer.writerow([gene.name,gene.uniprot_ID,gene.KEGG_ID,gene.MSI, gene.HRD,gene.FANC])