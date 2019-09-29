def get_expression(file_name):
    fh = open(file_name)
    global dict, col1, col2, col3, col4
    col1, col2, col3 = fh.readline().strip().split('\t')
    col4 = 'Exp'
    dict = {}
    lines = fh.readlines()
    for i in lines:
        x=i.rstrip().split('\t')
        if ((int(x[2])==0) or (int(x[1])==0)):
            exp = 's'
        else:
            if float(int(x[1])/int(x[2])) >= 1.5:
                exp = '+'
            elif float(int(x[2])/int(x[1]) )>= 1.5:
                exp = '-'
            else:
                exp = '.'
    # You insert the concatenated string as the dictionary value
            dict[x[0]]=x[1]+':'+x[2]+':'+exp
    return dict, col1, col2, col3, col4
def common_genes(dict1,dict2):
    global union_set
    set1 = set(list(dict1.keys()))
    set2 = set(list(dict2.keys()))
    union_set=set.union(set1,set2)
    return(union_set)
def unique_genes(dict1,dict2):
    global uniq_g1,uniq_g2
    set1 = set(list(dict1.keys()))
    set2 = set(list(dict2.keys()))
    uniq_g1 = set1 - set2
    uniq_g2 = set2 - (set1)
    print(uniq_g2)
    print(uniq_g1)
    return(uniq_g2,uniq_g1)
get_expression('cancer_lung_expression.txt')
dict1=dict
get_expression('cancer_prostate_expression.txt')
dict2=dict
common_genes(dict1,dict2)

unique_genes(dict1,dict2)

def compare():
