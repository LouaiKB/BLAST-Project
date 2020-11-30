"""

AUTHOR: Elise Fabre

"""
from ftplib import FTP
from sh import gunzip
from tkinter import * 
import pandas as pd
import string
import os.path

#function to search a specific oragnism name 
def search_organism(): 
    
    #open summary as a dataframe and sort it in order to get organism name and infraspecific_name with index summary_accesion number
    df = pd.read_table('assembly_summary_refseq.txt',header = 1,  sep = '\t')
    new_df = df[['organism_name','infraspecific_name']]
    
    #replace Nan values by empty boxes and merge 2 columns in one
    new_df['infraspecific_name'] = new_df['infraspecific_name'].fillna('')
    new_df["organism"]=new_df['organism_name']+ ' ' + new_df['infraspecific_name']
    new_df1 = new_df[['organism']]
    
    #create a new dataframe with lines that contain the chosen organism name 
    chosen_organism = entry_1.get()
    new_df2 = new_df1.loc[new_df1['organism'].str.contains(chosen_organism, regex=False, case=False, na=False)]
    
    #print the final dataframe in a label in the window with only 10 Proteoms
    label_list.config(text=new_df2.head(10))
    
    
#function to get proteom path from ftp server NCBI data base 
def proteom_path() :
    
    df = pd.read_table('assembly_summary_refseq.txt',header = 1,  sep = '\t')
    path_df = df[['organism_name','infraspecific_name', 'ftp_path']]
    
    #replace Nan values by empty boxes 
    path_df['infraspecific_name'] = path_df['infraspecific_name'].fillna('')
    path_df["organism"]=path_df['organism_name']+ ' ' + path_df['infraspecific_name']
    path_df2 = path_df[['organism', 'ftp_path']]

    #link chosen_accesion = entry_2 to the right ftp path 
    chosen_proteom = entry_2.get()
    index_nb= int(chosen_proteom)
    ftp_path= path_df2.loc[index_nb,"ftp_path"]
    path_part1 = ftp_path.rsplit("genomes/",1)[-1]+"/"
    path_part2 = ftp_path.rsplit('/')[-1]
    path_part3 = "_translated_cds.faa.gz"
    
    #ftp.path to download a specific file 
    #example :"all/GCF/000/002/285/GCF_000002285.3_CanFam3.1/GCF_000002285.3_CanFam3.1_translated_cds.faa.gz"
    full_path = path_part1+path_part2+path_part3
    
    #link chosen_accession to right organism and remove all punctuation
    path_organism = path_df2.loc[index_nb,'organism']
    path_organism = path_organism.translate(str.maketrans("","", string.punctuation))
    

    return (path_organism, full_path)
   

  

#function to download the chosen proteom 
def download_proteom():
    
    prot_path=proteom_path()
    
    #connection to ftp server NCBI
    ftp = FTP("ftp.ncbi.nih.gov")
    ftp.login()
    
    #Change directory
    directory = ("genomes/")
    ftp.cwd(directory)
    
    filename = prot_path[1]
    organism_path=prot_path[0]
    
    #download file and rename the downoalded file with oganism name
    with open(organism_path+".faa.gz", 'wb') as f:
           ftp.retrbinary('RETR ' + filename, f.write)
    
    #gunzip the downloaded file        
    return gunzip(organism_path+".faa.gz") 

# this function is created in order to check wether the file is downloaded or not 
def check_file():
    
    prot_path=proteom_path()
    organism_path=prot_path[0]
    
    #Check if a file has been downloaded to the working directory
    if os.path.isfile(organism_path+".faa"):
        label_3.configure(text="Yes, the proteome has been uploaded successfully !")
    else:
        label_3.configure(text="No, an error has occurreder")

              
#create a window to download a proteom
window = Tk()
window.title("Download proteom")
window.geometry("1000x500")

#first label 
label_1= Label(window, text = "Please enter an organism name : ")

#label to show proteom of chosen organism (max 10)
label_list = Label(window, justify=tkinter.LEFT,width=60, bg='SkyBlue2')

#second label 
label_2 = Label(window, text="Please enter the organism code :")

#Message label to confirm the download
label_3 = Label(window,justify=tkinter.LEFT)

#first entry : fild to write the searched organism name 
org_name = StringVar()
entry_1 = Entry(window,textvariable=org_name, width=10, bg='SkyBlue2')

#fild to write the accession number of the choosen organism
accession = StringVar()
entry_2 = Entry(window,textvariable=accession, width=10, bg='SkyBlue2')

#first button to display the list that contains the choosen organism in the label
button_1 = Button(window, text="Search proteoms", command = search_organism)

#second button to download the proteom 
button_2 = Button(window, text="Download a proteom", command=download_proteom)

#third button to check if the file is in the working directory
button_3 = Button(window, text="Is the proteome downloaded?", command=check_file)


label_1.grid(row=1,column=1)
label_list.grid(row=2, column=3)
label_2.grid(row=4, column=1)
label_3.grid(row=6, column=3)

entry_1.grid(row=1, column=2)
entry_2.grid(row=4, column=2)

button_1.grid(row=1, column=3)
button_2.grid(row=4, column=3)
button_3.grid(row=5, column=3)


window.mainloop()
window.destroy()