#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 14:21:37 2022

@author: wenrchen
"""

#match uniprot ptm annotation with proteoforms

import pandas as pd
import sys
import warnings
warnings.filterwarnings("ignore")

def extract_gene(tmp):
    # gene_name=tmp.split("|")[-2] // original
    gene_name = tmp.split('|')[-1].split('_')[0]
    return gene_name

def checkInt(tmp):
    try:
        int(tmp)
        return True
    except ValueError:
        return False

def get_ptm_range(tmp,first):
    #remove the letter before "." and the letter after "."

    #remove the letter between "[" and  "]"(include"[" and "]")
    while(tmp.find("[")!=-1):
        
        if(tmp.find("[Acetyl]")!=-1):
            flag=tmp.find("[Acetyl]")
            tmp=tmp[:flag-3]+tmp[flag-2]+tmp[flag+len("[Acetyl]"):]
    
        while(tmp.find("(C)[Carbamidomethylation]")!=-1):
            flag=tmp.find("(C)[Carbamidomethylation]")
            tmp=tmp[:flag]+tmp[flag+1]+tmp[flag+len("(C)[Carbamidomethylation]"):]
        
        flag1=tmp.find("[")
        flag2=tmp.find("]")
        tmp=tmp[:flag1]+tmp[flag2+1:]
        
    flag=tmp.find(".")
    tmp=tmp[flag+1:]
    flag=tmp.find(".")
    tmp=tmp[:flag]
    
    # while(tmp.find("(C)")!=-1):
    #     flag=tmp.find("(C)")
    #     tmp=tmp[:flag]+tmp[flag+1]+tmp[flag+3:]
        
    
    flag1=tmp.find("(")
    flag2=tmp.find(")")

    return flag1+first, flag2+first-2

def mass_filter(uniprot_evidence, df_ptm):
    entries = uniprot_evidence['Protein accession'].apply(lambda x: x.split('|')[1])
    uniprot_idx = uniprot_evidence.index.values
    mass_v_list = []
    idx_ext = []
    for i in range(len(entries)):
        mass_df = df_ptm.loc[df_ptm['Entry'] == entries[i], :]
        if len(mass_df):
            mass_v = mass_df['Mass'].values[0]
            mass_v_list.append(mass_v)
            idx_ext.append(uniprot_idx[i])

    uniprot_evidence = uniprot_evidence.loc[idx_ext, :]
    uniprot_evidence.insert(21,"Uniprot mass", mass_v_list)
    # filter
    idx_list = uniprot_evidence.index.values
    idx_ext = []
    tol = 1.2
    for i in range(len(uniprot_evidence)):
        diff_mass = uniprot_evidence['Mass shift'].iloc[i] - uniprot_evidence['Uniprot mass'].iloc[i]
        if abs(diff_mass) < tol:
            idx_ext.append(idx_list[i])

    uniprot_evidence_ext = uniprot_evidence.loc[idx_ext, :]
    return uniprot_evidence_ext

df=pd.read_csv(sys.argv[1],sep='\t')
#df=pd.read_excel(sys.argv[1],sheet_name=0)

df_ptm=pd.read_csv(sys.argv[2],sep='\t')
output=sys.argv[3]

df_ptm=df_ptm.dropna()
print(df_ptm.shape[0])

ptm_dict={}
for i in range(df_ptm.shape[0]):
    gene=df_ptm.iloc[i]['Gene names  (primary )'].split(";")[0]
    ptms=df_ptm.iloc[i]['Modified residue'].split("MOD_RES ")[1:]
    
    if(gene not in ptm_dict.keys()):
        ptm_dict[gene]=ptms
    else:
        for ptm in ptms:
            ptm_dict[gene].append(ptm)

uniprot_evidence=pd.DataFrame(columns=df.columns)
uniprot_ptms=[]
matched_index=[]
#potential_C57=[]

for i in range(df.shape[0]):
    protein_id=df.iloc[i]['Protein accession']
    gene_name=extract_gene(protein_id)
    # print(gene_name)
    # print(ptm_dict.keys())
    
    if(gene_name in ptm_dict.keys()):
        first_residue=df.iloc[i]['First residue']
        last_redisue=df.iloc[i]['Last residue']
        
        ptms=ptm_dict[gene_name]
        for ptm in ptms:
            pos=ptm.split(";")[0]
            # print(pos)
            # print(ptm)
            if(checkInt(pos)):
                pos=int(pos)
                if(pos>=first_residue and pos<=last_redisue):
                    seq=df.iloc[i]['Proteoform']
                    ptm_first,ptm_last=get_ptm_range(seq, first_residue)
                    #ptm_first,ptm_last=first_residue,first_residue
                    if(pos>=ptm_first and pos<=ptm_last):
                    #     if(abs(df.iloc[i]['Mass shift']+57)<=1):
                    #         potential_C57.append(i)
                        #if(ptm.find("N-acetyl")!=-1):
                            matched_index.append(i)
                            uniprot_evidence=uniprot_evidence.append(df.loc[i],ignore_index=True)
                            uniprot_ptms.append(ptm)

print(len(list(set(matched_index))))
# print(len(list(set(potential_C57))))

uniprot_evidence.insert(20,"Uniprot evidence", uniprot_ptms)
# mass filter
uniprot_evidence_ext = mass_filter(uniprot_evidence, df_ptm)
uniprot_evidence_ext.to_csv(output,sep='\t',index=None)
        
        
    
