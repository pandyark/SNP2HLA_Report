import ast
import json
from ast import literal_eval

import length as length
import pandas as pd
from HLA_Dictionary import HLA_dict


data = pd.read_csv("/Users/bowcock_lab/Downloads/SNP2HLA_package_v1.0.3/SNP2HLA/Psoriasis_EUR_IMPUTED.bgl.phased", sep="\t", header= None )
data = data.drop(data.columns[0], axis = 1)

# Will print the first column of the pandas dataframe

#print(data.loc[0:,1])
data_length = (len(data))

#convert hla_dict to list
hla_values_list=[]
for key, value in HLA_dict.items():
    temp = value
    for item in temp:
        hla_values_list.append(item)
hla_values_list.append("pedigree")
hla_values_list.append("id")
#data[data['1'].isin()]

selected_hla_data = data[data[1].isin(hla_values_list)]
#selected_hla_data = data[data[1].]
selected_hla_data.to_csv("/Users/bowcock_lab/Desktop/hla_test.csv")

sample_ids = data.loc[0,2:]
mylist = list(dict.fromkeys(sample_ids))
final_values_for_report = []
columns = ["Sample_ID", "HLA_A_1", "HLA_A_2", "HLA_C_1", "HLA_C_2", "HLA_B_1", "HLA_B_2"]
for i in mylist:
    new_hla_values =[]
    columns_to_process = []
    row_vals = list(selected_hla_data.loc[1,:])
    columns_to_process.append(1)

    for j in range(0, len(row_vals)):
        if i==row_vals[j]:
            columns_to_process.append(j+1)
    data_to_work =  selected_hla_data
    # print(columns_to_process)
    dataset = data_to_work[data_to_work.columns.intersection(columns_to_process)]
    newdataset = dataset.drop(dataset.index[1])
    final = newdataset[(newdataset[columns_to_process[1]] == "P") | (newdataset[columns_to_process[2]] == "P") | (newdataset[columns_to_process[1]] == i)].transpose()
    #hla_report_final.append(final.iloc[0])
    new_hla_values.append(i)
    new_hla_values.extend(list(final.iloc[0]))

    new_hla_values.remove("pedigree")

    if len(new_hla_values) < 7:

        new_hla_values.extend('NA' for i in range (7- len(new_hla_values)))

        if not "HLA_A" in new_hla_values[1]:
            new_hla_values.insert(1,new_hla_values.pop())

        if not "HLA_A" in new_hla_values[2]:
            new_hla_values.insert(2, new_hla_values.pop())

        if not "HLA_C" in new_hla_values[3]:
            new_hla_values.insert(3, new_hla_values.pop())

        if not "HLA_C" in new_hla_values[4]:
            new_hla_values.insert(4, new_hla_values.pop())

        final_values_for_report.append(new_hla_values)

    if not len(new_hla_values) > 7:
        final_values_for_report.append(new_hla_values)



#print(len(final_values_for_report))
hla_report_final = pd.DataFrame(final_values_for_report, columns=["Sample_ID", "HLA_A_1", "HLA_A_2", "HLA_C_1", "HLA_C_2", "HLA_B_1", "HLA_B_2"])

#print(hla_report_final)
hla_report_final.to_csv("/Users/bowcock_lab/Desktop/hla_final_report.csv")






