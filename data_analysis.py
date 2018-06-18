#read in dna from text file

def read_dna(dna_file):
  dna_data=""
  count=0
  with open(dna_file, "r") as file:
    for line in file:
      if (line[0:1] !='>') and count < 7000:
        dna_data+=(line.strip('\n')).upper()
        count+=1
        #dna_data.append(line.strip('\n'))
     
  return dna_data

#separates dna into strands of specified length

def generate_samples(dna_data):
  sample_list=[];
  for i in range(0,len(dna_data),300):
    sample_list.append(dna_data[i:i+300])  
     
  return sample_list 

#Separates DNA strand into codons of length num

def dna_codons(strand_sample,num):
    strand_sample_sep=[strand_sample[i:i+num] for i in range(0,len(strand_sample)-num+1)]
    
    return strand_sample_sep


def nt_dict_array(codon_list,num):
#given the codon length, counts frequency for each permutation in fragment
#returns frequency of each permutation in the fragment as array
  import numpy as np  
  from collections import Counter
  nt=['A','C','G','T']
  r1=range(0,4)  
  if num ==1:
    perm=['A','C','G','T']
  elif num == 2:
    perm=[nt[i]+nt[j] for i in r1 for j in r1]  
  elif num == 3:
    perm=[nt[i]+nt[j]+nt[k] for i in r1 for j in r1 for k in r1]  
  elif num == 4:
    perm=[nt[i]+nt[j]+nt[k]+nt[l] for i in r1 for j in r1 for k in r1 for l in r1]  
  sample_array=np.zeros(4**num)
  lst=Counter(codon_list)
  for i in range(0,4**num):
    sample_array[i]=lst[perm[i]]
  return sample_array
  
def Xmatrix_generator(dna_data,num,filename):
#generates matrix (mxn array) with frequency of codons (m = 4**num total)... 
#for each fragment (n total)
  import csv
  import numpy as np
  Xmatrix=np.zeros((4**num,len(dna_data)))
  for i in range(0,len(dna_data)):
    Xmatrix[:,i]= nt_dict_array(dna_data[i],num)  
  myFile = open(filename, 'w')
  with myFile:
    writer = csv.writer(myFile)
    writer.writerows(Xmatrix)  
  return Xmatrix  

def Seq_Analysis(readfile,num,outfile):
#import in the dna file and separates it into fragment lengths of 300
#for each fragment, it creates codons of length num (num = 3 in real dna)
  import csv
  import numpy as np
  x=read_dna(readfile)
  strand_samples=generate_samples(x)
  dna_data=[dna_codons(strand_samples[i],num) for i in range(0,len(strand_samples))]
  Xmatrix_generator(dna_data,num,outfile)    

#Reads in file and outputs frequencies

import numpy as np
import csv

#Genomic Data and Codon length
file_name='Pseudomonas_aeruginosa_genomic.txt'
num=3

Seq_Analysis(file_name,num,'Xmatrix_num3pa.csv')

#reads in frequency matrices for fragments in each case
reader = csv.reader(open("Xmatrix_num3pa.csv", "r"), delimiter=",")
x = list(reader)
result = np.array(x).astype("float")


#Data to be used later
x=read_dna('Pseudomonas_aeruginosa_genomic.txt') 
strand_samples=generate_samples(x)

from sklearn.decomposition import PCA as sklearnPCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

#Initial PCA for 2d plot
sklearn_pca = sklearnPCA(n_components=2)

#Initialize PCA for 3d plot
sklearn_pcahigh = sklearnPCA(n_components=3)

#PCA transformation for each codon length in 2D and 3D
sklearn_transf = sklearn_pca.fit_transform(result.T)
sklearn_transfhigh = sklearn_pcahigh.fit_transform(result.T)


#Find Centroids (Verts) of Data After PCA
kmeans = KMeans(n_clusters=7,random_state=0).fit(sklearn_transfhigh)
verts=kmeans.cluster_centers_

#Determine End Codon Count in Each Cluster
from collections import Counter

def end_codon_count(end_codon,codon_list,cluster_pred):
  from operator import itemgetter
  
  #creates a list of tuple of end codon frequency and cluster number
  lst=[(codon_list[i].count(end_codon),cluster_pred[i]) for i in range(0,len(codon_list))]
  
  #sorts lst by cluster number so that we can find total end codon count
  lst.sort(key=itemgetter(1))
  
  #makes into array to perform operations on columns and rows
  lstarray=np.array(lst)
  clind=[np.where(lstarray[:,1]==j)[0][1] for j in range(7)]
  clind.append(len(lst)-1)
  
  #calculates the mean instances of end codon
  meanend_codon_count= [np.sum(lstarray[clind[i]:clind[i+1],0])/(clind[i+1]-clind[i]) for i in range(len(clind)-1)] 
  
  #calculates the total instances of end codon
  end_codon_count=[np.sum(lstarray[clind[i]:clind[i+1],0]) for i in range(len(clind)-1)] 
  
  return np.array(meanend_codon_count)
  
#calculate the mutual information for determine the 

r1=range(0,4)
nt=['A','C','G','T']
x=read_dna(file_name)
strand_samples=generate_samples(x)
dna_data=[dna_codons(strand_samples[i],3) for i in range(0,len(strand_samples))]
  
  
perm3=[nt[i]+nt[j]+nt[k] for i in r1 for j in r1 for k in r1]  
ldna=len(dna_data[0])

cluster_mi=[]

#calculates the mutual information for each cluster
for k in range(len(dna_data)):
  lst=Counter(dna_data[k])
  sumList1=[];


#Finds frequencies of each nucleotide in fragment p1i, p2j, p3k  
  for j in range(4):
    sumList1.append(sum([lst[permlist1[j][i]] for i in range(16)]))

  d1=dict(zip(nt,sumList1))

  sumList2=[];
  for j in range(4):
    sumList2.append(sum([lst[permlist2[j][i]] for i in range(16)]))
  d2=dict(zip(nt,sumList2))

  sumList3=[];
  for j in range(4):
    sumList3.append(sum([lst[permlist3[j][i]] for i in range(16)]))
  d3=dict(zip(nt,sumList3))    
    
  mi_sum=0
  
  for i in range(64):
   if lst[perm3[i]]>0:   
     mi_sum+=lst[perm3[i]]*np.log(ldna**2*lst[perm3[i]]/(d1[perm3[i][0]]*d2[perm3[i][1]]*d3[perm3[i][2]]))/ldna
  
  cluster_mi.append((mi_sum,y_pred[k]))

from operator import itemgetter

cluster_mi.sort(key=itemgetter(1))
clustmi_array=np.array(cluster_mi)
clusterind=[np.where(clustmi_array[:,1]==j)[0][0] for j in range(7)]
clusterind.append(len(cluster_mi)-1)
mean_mi=[np.sum(clustmi_array[clusterind[i]:clusterind[i+1],0])/(clusterind[i+1]-clusterind[i]) for i in range(len(clusterind)-1)] 

