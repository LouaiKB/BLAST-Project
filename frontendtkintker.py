from clustersGenerator import ClustersGenerator
from tkinter import *
from tkinter import filedialog,messagebox
import os 

list_files=[]


def browse_file():
    file =filedialog.askopenfilenames(title="Select at least 2 files",
                                                filetypes=[( "fasta file" , ".fasta" ), 
                                                          ( "faa" , ".faa" )])
    global list_files 
    for i in file:
        list_files.append(i)
    label_2.configure(text="\n".join(list_files))


def run_blast():
    lablel_3.configure(text='In process .....')
    global list_files
    new_list_files = []
    for file in list_files:
        fastaName = file.split("/")[-1].replace(' ', '')
        new_list_files.append(fastaName)
    for i in range(len(new_list_files)) :    
         os.rename(list_files[i],new_list_files[i])
    blast_object=ClustersGenerator(new_list_files)
    blast_object.clusterization()
    


#create awindow to run blast 
window = Tk()
window.title("Make clusters of orthologs proteins")
window.geometry("1100x600")

#first label 
label_1= Label(window, text = "Please choose at leat 2 proteoms : ")
#label to show different file path selected
label_2=Label(window,text="opened files:", width=100)

#first button to browse files 
button_1=Button(window, text="browse files",command=browse_file)
#Second button to launch blast process 
button_2=Button(window, text="Blast and clusters generator", command=run_blast)
lablel_3 = Label(window, text='')

#button to close the window
button_quit=Button(window, text="Quit", command=window.destroy)

label_1.grid(row=1, column=1)
label_2.grid(row=2, column=2)
lablel_3.grid(row=5, column=2)

button_1.grid(row=1, column=2)
button_2.grid(row=6, column=2)
button_quit.grid(row=15, column=2)



window.mainloop()
# window.destroy()