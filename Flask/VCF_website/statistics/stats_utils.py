import ast

def decompress(gt_arr):
    """ Data for geneotyype array is stored as compressed data using similar technique to Run-length encoding (RLE)
    the compressed data includes list of lists [['2','a']]. This corresponds to [[0,0],[0,0]] : two samples that are
    homozygous alternative. The positions of the sample are also saved.  To decompress we use the value that corespond to 
    the compressed value and multiply with the integer
    """
    freq_data = []
    gt_data = []
    decomp_dict = {'a':'[0, 0]','b':'[0, 1]','c':'[1, 0]','d':'[1,1]'}
    for item in gt_arr:
        gt_arr_data = ast.literal_eval(item['genotypes'])
        freq_data.append(ast.literal_eval(item['geno_freq']))

        snp_data = ""
        for sample in gt_arr_data:
            rep = int(sample[0])
            val = decomp_dict[sample[1]]+','
            val = rep*val
            snp_data += val

        gt_data.append(list(ast.literal_eval(snp_data[:-1])))

    return gt_data, freq_data