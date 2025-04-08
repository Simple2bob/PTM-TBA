#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 10:48:24 2023

@author: wenrchen
"""
#match dbPTM ptm annotation with proteoforms

import pandas as pd
import sys
import warnings
warnings.filterwarnings("ignore")

def extract_gene(tmp):
    gene_name=tmp.split("|")[-2] ##-2// original
    # gene_name = tmp.split('|')[-1].split('_')[0]
    return gene_name

def checkInt(tmp):
    try:
        int(tmp)
        return True
    except ValueError:
        return False
    
fix_mod_form1="(C)[Carbamidomethylation]"
fix_mod_form2="C[Carbamidomethylation]"

def get_ptm_range(tmp,first,n_term_flag,variable_flag):
    #remove the letter before "." and the letter after "."

    #remove the letter between "[" and  "]"(include"[" and "]")
    
    if(n_term_flag==1):
        flag=tmp.find("[42.0106]")
        if(tmp[flag-1]==")"):
            tmp=tmp[:flag-3]+tmp[flag-2]+tmp[flag+9:]
        else:
            tmp=tmp[:flag]+tmp[flag+9:]
    
    while(tmp.find("[")!=-1):
        
        
        while(tmp.find(fix_mod_form1)!=-1):
            flag=tmp.find(fix_mod_form1)
            tmp=tmp[:flag]+tmp[flag+1]+tmp[flag+len(fix_mod_form1):]
        
        while(tmp.find(fix_mod_form2)!=-1):
            flag=tmp.find(fix_mod_form2)
            tmp=tmp[:flag+1]+tmp[flag+len(fix_mod_form2):]
        
        flag1=tmp.find("[")
        flag2=tmp.find("]")
        tmp=tmp[:flag1]+tmp[flag2+1:]
        
        if(tmp.find("[")==-1):
            
            if(variable_flag==1):
                if(tmp[flag1-1]!=")"):
                    tmp=tmp[:flag1-1]+"("+tmp[flag1-1]+")"+tmp[flag1:]

    flag=tmp.find(".")
    tmp=tmp[flag+1:]
    flag=tmp.find(".")
    tmp=tmp[:flag]
    
    flag1=tmp.find("(")
    flag2=tmp.find(")")

    return flag1+first,flag2+first-2


def mass_filter(dbptm_evidence):
    # common ptms
    ptm_mass_list = {'Acetylation': 42.010565,
                     'Phosphorylation': 79.966331,
                     'Methylation': 14.01565,
                     'Oxidation': 15.994915}
    dbptm_evidence['dbPTM evidence'] = dbptm_evidence['dbPTM evidence'].astype(str)
    dbptm_evidence['PTM name'] = dbptm_evidence['dbPTM evidence'].str.extract(r"\(([^,]+),")
    ptm_type = dbptm_evidence['PTM name'].unique()[0]
    ptm_type = ptm_type.strip(" '\'")  # remove quotes and spaces
    if ptm_type not in ptm_mass_list:
        print(f"{ptm_type} is not available, only annotations with Acetylation, Phosphorylation, Methylation, Oxidation provided!")
        return
    # filter
    tol = 1.2
    dbptm_evidence_ext = dbptm_evidence[abs(dbptm_evidence['Mass shift'].values - ptm_mass_list[ptm_type]) < tol]
    dbptm_evidence_ext = dbptm_evidence_ext.drop(columns="PTM name", axis=1)
    return dbptm_evidence_ext

#df=pd.read_csv(sys.argv[1],sep='\t') 
#df=pd.read_excel(sys.argv[1],sheet_name=0)
dfname=sys.argv[1]        
if(dfname.find("xlsx")==-1):
   df=pd.read_csv(dfname,sep='\t') ##proteoforms with PTM df
else:
    df=pd.read_excel(dfname,sheet_name=0)

df_ptm=pd.read_csv(sys.argv[2],sep='\t',header=None)
output=sys.argv[3]
n_term_mode=int(sys.argv[4])  #0 means"no n-term"

# gene_uniprot_df=pd.read_csv(sys.argv[4],sep='\t')
# gene_uniprot_df = gene_uniprot_df.dropna(subset=["Gene names  (primary )", "Entry"])
# gene_uniprot_map={}
# for i in range(gene_uniprot_df.shape[0]):
#     # gene_name=gene_uniprot_df.iloc[i]["From"] // original
#     # gene_name=gene_uniprot_df.iloc[i]["Gene names  (primary )"]
#     gene_name=gene_uniprot_df.iloc[i]["Entry"]
#
#     gene_uniprot_map[gene_name]=gene_uniprot_df.iloc[i]["Entry"]

df_ptm=df_ptm.dropna()
print(df_ptm.shape[0])

ptm_dict={}
for i in range(df_ptm.shape[0]):
    gene=df_ptm.iloc[i][1]
    ptm_pos=df_ptm.iloc[i][2]
    ptm=df_ptm.iloc[i][3]
    ptm_pubid=df_ptm.iloc[i][4]
    
    
    
    if(gene not in ptm_dict.keys()):
        ptm_dict[gene]=[(ptm,ptm_pos,ptm_pubid)]
    else:
        ptm_dict[gene].append((ptm,ptm_pos,ptm_pubid))

dbptm_evidence=pd.DataFrame(columns=df.columns)
dbptm_ptms=[]
matched_index=[]

for i in range(df.shape[0]):
    protein_id=df.iloc[i]['Protein accession']
    gene_name=extract_gene(protein_id)
    # uniprot_id=gene_uniprot_map[gene_name]
    
    # if(uniprot_id in ptm_dict.keys()):
    if (gene_name in ptm_dict.keys()):
        first_residue=df.iloc[i]['First residue']
        last_redisue=df.iloc[i]['Last residue']
        
        # ptms=ptm_dict[uniprot_id]
        ptms=ptm_dict[gene_name]
        for ptm in ptms:
            pos=ptm[1]
            if(checkInt(pos)):
                pos=int(pos)
                if(pos>=first_residue and pos<=last_redisue):
                    seq=df.iloc[i]['Proteoform']
                    
                    n_term_flag=0
                    if(df.iloc[i]['Protein N-terminal form'].find("ACETYLATION")!=-1):
                        n_term_flag=1
                    variable_flag=0
                    if(df.iloc[i]['#variable PTMs']==1):
                        variable_flag=1
                    
                    if(n_term_mode==1):
                        ptm_first=first_residue
                        if(pos==ptm_first):
                            matched_index.append(i)
                            dbptm_evidence=dbptm_evidence.append(df.loc[i],ignore_index=True)
                            dbptm_ptms.append(ptm)
                    else:
                        ptm_first,ptm_last=get_ptm_range(seq, first_residue,n_term_flag,variable_flag)
                        if(pos>=ptm_first and pos<=ptm_last):
                            matched_index.append(i)
                            dbptm_evidence=dbptm_evidence.append(df.loc[i],ignore_index=True)
                            dbptm_ptms.append(ptm)

print(len(list(set(matched_index))))

dbptm_evidence.insert(20,"dbPTM evidence", dbptm_ptms)
dbptm_evidence_ext = mass_filter(dbptm_evidence)
dbptm_evidence_ext.to_csv(output,sep='\t',index=None)