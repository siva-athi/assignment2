from Bio import SeqIO

file_handler=open('cancer_lung_expression.txt')
file=file_handler.readlines()[1:]
for i in file:
    a=i.rstrip().split('\t')
    gene_name=a[0]
    num_l= a[1]
    num_c= a[2]
    if float(int(num_l)/int(num_c)>=1.5):
        print(str(gene_name)+'\t'+str(num_l)+':'+str(num_c)+':'+'+')
    elif float(int(num_c) / int(num_l) >= 1.5):
        print(str(gene_name) + '\t' + str(num_l) + ':' + str(num_c)+':'+'-')
    elif float(int(num_l) / int(num_c)== 1):
        print(str(gene_name) + '\t' + str(num_l) + ':' + str(num_c)+':'+'.')