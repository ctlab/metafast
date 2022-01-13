import numpy as np
import pandas as pd
from scipy.spatial import distance
from seaborn import clustermap

Air = ["taxonomic_profile_4.txt" ,"taxonomic_profile_7.txt" ,"taxonomic_profile_8.txt" ,"taxonomic_profile_9.txt" ,"taxonomic_profile_10.txt" ,"taxonomic_profile_11.txt" ,"taxonomic_profile_12.txt" ,"taxonomic_profile_23.txt" ,"taxonomic_profile_26.txt" ,"taxonomic_profile_27.txt"]  
Skin = ["taxonomic_profile_1.txt" ,"taxonomic_profile_13.txt" ,"taxonomic_profile_14.txt" ,"taxonomic_profile_15.txt" ,"taxonomic_profile_16.txt" ,"taxonomic_profile_17.txt" ,"taxonomic_profile_18.txt" ,"taxonomic_profile_19.txt" ,"taxonomic_profile_20.txt" ,"taxonomic_profile_28.txt"]  
Oral = ["taxonomic_profile_6.txt" ,"taxonomic_profile_7.txt" ,"taxonomic_profile_8.txt" ,"taxonomic_profile_13.txt" ,"taxonomic_profile_14.txt" ,"taxonomic_profile_15.txt" ,"taxonomic_profile_16.txt" ,"taxonomic_profile_17.txt" ,"taxonomic_profile_18.txt" ,"taxonomic_profile_19.txt"]  
Gas = ["taxonomic_profile_0.txt" ,"taxonomic_profile_1.txt" ,"taxonomic_profile_2.txt" ,"taxonomic_profile_3.txt" ,"taxonomic_profile_4.txt" ,"taxonomic_profile_5.txt" ,"taxonomic_profile_9.txt" ,"taxonomic_profile_10.txt" ,"taxonomic_profile_11.txt" ,"taxonomic_profile_12.txt"]  
Uro = ["taxonomic_profile_0.txt" ,"taxonomic_profile_2.txt" ,"taxonomic_profile_3.txt" ,"taxonomic_profile_5.txt" ,"taxonomic_profile_6.txt" ,"taxonomic_profile_21.txt" ,"taxonomic_profile_22.txt" ,"taxonomic_profile_24.txt" ,"taxonomic_profile_25.txt"]  
os.chdir("./Air/short_read")
C = pd.DataFrame()
for i in range(len(Air)):
    tmp = pd.read_csv(Air[i], header=3, delimiter="\t")
    C = C.append(tmp[tmp['RANK']=='species']['PERCENTAGE'], ignore_index=True)
os.chdir("./Skin/short_read")
U1 = pd.DataFrame()
for i in range(len(Skin)):
    tmp = pd.read_csv(Skin[i], header=3, delimiter="\t")
    U1 = U1.append(tmp[tmp['RANK']=='species']['PERCENTAGE'], ignore_index=True)
C = C.append(U1, ignore_index=True)
os.chdir("./Oral/short_read")
U1 = pd.DataFrame()
for i in range(len(Oral)):
    tmp = pd.read_csv(Oral[i], header=3, delimiter="\t")
    U1 = U1.append(tmp[tmp['RANK']=='species']['PERCENTAGE'], ignore_index=True)
C = C.append(U1, ignore_index=True)
os.chdir("./Gastro/short_read")
U1 = pd.DataFrame()
for i in range(len(Gas)):
    tmp = pd.read_csv(Gas[i], header=3, delimiter="\t")
    U1 = U1.append(tmp[tmp['RANK']=='species']['PERCENTAGE'], ignore_index=True)
C = C.append(U1, ignore_index=True)
os.chdir("./Uro/short_read")
U1 = pd.DataFrame()
for i in range(len(Uro)):
    tmp = pd.read_csv(Uro[i], header=3, delimiter="\t")
    U1 = U1.append(tmp[tmp['RANK']=='species']['PERCENTAGE'], ignore_index=True)
C = C.append(U1, ignore_index=True)
X = C.fillna(0).values

dist_matrix = np.empty((X.shape[0],X.shape[0]),dtype=float) 
for i in range(X.shape[0]):
    for j in range(X.shape[0]):
        dist_matrix[i,j] = distance.braycurtis(X[i],X[j])
      
rcols= ["Airways_4" ,"Airways_7" ,"Airways_8" ,"Airways_9" ,"Airways_10" ,"Airways_11" , "Airways_12" ,"Airways_23" ,"Airways_26" ,"Airways_27","Skin_1" ,"Skin_13" ,"Skin_14" ,"Skin_15" ,"Skin_16" ,"Skin_17" ,"Skin_18" ,"Skin_19" ,"Skin_20" ,"Skin_28","Oral_6" ,"Oral_7" ,"Oral_8" ,"Oral_13" ,"Oral_14" ,"Oral_15" ,"Oral_16" ,"Oral_17" ,"Oral_18" ,"Oral_19", "Gastro_0" ,"Gastro_1" ,"Gastro_2" ,"Gastro_3" ,"Gastro_4" ,"Gastro_5" ,"Gastro_9" ,"Gastro_10" ,"Gastro_11" ,"Gastro_12" , "Uro_0" ,"Uro_2" ,"Uro_3" ,"Uro_5" ,"Uro_6" ,"Uro_21" ,"Uro_22" ,"Uro_24" ,"Uro_25"]  
dist_df = pd.DataFrame(result,index=rcols,columns=rcols)
clustermap(dist_matrix, method ='ward', figsize=(16,16))

