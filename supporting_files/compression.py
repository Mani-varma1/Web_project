import ast
import pandas as pd


def compress(array):
    diict = {'[0, 0]':'a','[0, 1]':'b','[1, 0]':'c','[1, 1]':'d'}
    """ Get the actual data from the array to wrangle"""
    arr = ast.literal_eval(array)
    
    # turn it one sample genotype array into  a string and compress it 
    temp = []
    for item in arr:
        value = diict[str(item)]
        temp.append(value)
        
    compressed = []
    prev = None
    count = 0
    for i,ele in enumerate(temp):
        if prev == None:
            count += 1
            prev = ele

        elif ele == prev and i != len(temp)-1:
            count += 1

        elif ele != prev and i != len(temp)-1:
            lst = [count,prev]
            compressed.append(lst)
            count = 1
            prev = ele

        else:
            if ele == prev and i == len(temp)-1:
                count+=1
                lst = [count,prev]
                compressed.append(lst)
                count = 0 
                prev = None
            elif ele != prev and i == len(temp)-1:
                lst = [count,prev]
                compressed.append(lst)
                count =1
                prev = ele
                lst = [count,prev]
                compressed.append(lst)
                prev = None
                count = 0
                
    return compressed

def read_and_compress(input):
    df = pd.read_csv(input,sep='\t')
    gt_df = df.loc[:,'genotype_array']
    final_lst = []
    for item in gt_df:
        gt_arr =  compress(item)
        final_lst.append(gt_arr)
    df['genotype_array']= final_lst
    return df 

df = read_and_compress("YRI_data.tsv")
df.to_csv('YRI_compressed.tsv', sep = '\t')

df = read_and_compress("MXL_data.tsv")
df.to_csv('MXL_compressed.tsv', sep = '\t')

df = read_and_compress("GBR_data.tsv")
df.to_csv('GBR_compressed.tsv', sep = '\t')

df = read_and_compress("JPT_data.tsv")
df.to_csv('JPT_compressed.tsv', sep = '\t')

df = read_and_compress("PJL_data.tsv")
df.to_csv('PJL_compressed.tsv', sep = '\t')


