from Bio import SeqIO
import re
import sys
import getopt
argv = sys.argv[1:]
gene = ''
lung = ''
pros = ''
try:
    opts, args = getopt.getopt(argv, "hf:l:p:", ["faa=","lung=","pros="])
except getopt.GetoptError:
    print('Please enter invalid option(s)!')
    print('Usage: lakshas_a2.py -f <faa> -l <lung>,-p <prostate>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print('this files test gene expression.')  # need more work
        sys.exit()
    elif opt in ("-f", "--faa"):
        gene = arg
    elif opt in ("-l", "--lung"):
        lung = arg
    elif opt in ("-p", "--pros"):
        pros = arg
def get_seq_fasta(filename):
    global mydict
    record_iterator = SeqIO.parse(filename,"fasta")
    mydict={}
    for i in record_iterator:
        mydict[i.id] = str(i.seq)
    return(mydict)
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
    global common_set,s
    set1 = set(list(dict1.keys()))
    set2 = set(list(dict2.keys()))
    common_set=set.intersection(set1,set2)
    #print(common_set)
    s=""
    for i in common_set:
        s=s+i+', '
    return(common_set)
def unique_genes(dict1,dict2):
    global uniq_g1,uniq_g2,v,u
    set1 = set(list(dict1.keys()))
    set2 = set(list(dict2.keys()))
    uniq_g1 = set1 - set2
    u=""
    uniq_g2 = set2 - set1
    v=""
    #print(uniq_g1)
    for i in uniq_g1:
        u=u+i+', '
    #print(uniq_g2)
    for i in uniq_g2:
        v=v+i+', '
    return(uniq_g2,uniq_g1)
def compare(dict1, dict2,common_set):
    for i in common_set:
        c=i
        c+= " Numl:Numc:Exp "
        c+='('+dict1[i]+')'
        a = dict1[i].split(":")
        c += " NumP:Numc:Exp "
        c+='('+dict2[i]+')'
        b = dict2[i].split(":")
        #print(b[2)
        if a[2]==b[2]:
            c+='congurnece'
        else:
            c+='disparcity'
        print(c)
    return ()
def motif_finder(mydict):
    global dict123,dic456
    dict123 = {}
    for k,v in mydict.items():

        comp=re.compile("(AAG[ATGC]{0,2}AA)")
        if k in z:
            c= (comp.findall(v))
            dic456={}
            for i in c:
                if i not in dic456.keys():
                    dic456[i] = 1
                else:
                    dic456[i]+= 1
            dict123[k]= dic456
    #print(dict123)
    return(dict123)


get_seq_fasta(gene)

print("1. Files: FASTA "+"("+gene+")"+ " Cancer Exp "+ "("+lung+","+pros+")\n")
get_expression(lung)
dict1=dict
print('2. Expression Data for cancer_lung_expression.txt')
print(col1+" "+col2+":"+col3+":"+col4)
for k,v in dict1.items():
    print(k+' '+v)
print("...\n")
get_expression(pros)
dict2=dict
print('3. Expression Data for cancer_prostate_expression.txt')
print(col1+" "+col2+":"+col3+":"+col4)
for k,v in dict2.items():
    print(k+' '+v)
print("...\n")
common_genes(dict1,dict2)
print("4. Expression Comparison (cancer_lung_expression.txt) vs (cancer_prostate_expression.txt)")
compare(dict1,dict2,common_set)
print("...\n")
print("5. Common Genes for cancer_lung_expression.txt and cancer_prostate_expression.txt")
print(s[:-2])
print("\n")
unique_genes(dict1,dict2)
print("6. Genes unique for cancer_lung_expression.txt")
print(u[:-2])
print("\n")
print("7. Genes unique for cancer_prostate_expression.txt")
print(v[:-2])
print("\n")
q1=input("Which set do you want to examine (5|6|7)")
print("\n")
print("8. Your answer is [{0}]".format(q1))
print("We will do motif finding for the following sequences:")
if q1=="5":
    z=common_set
    print(s[:-2])
    motif_finder(mydict)
    for k,v in dict123.items():
        print(k)
        print("Motif\tFrequency")
        for g,h in v.items():
            print(g+"\t"+str(h))
        print("...\n")
elif q1=="6":
    z=uniq_g1
    print(u[:-2])
    motif_finder(mydict)
    for k,v in dict123.items():
        print(k)
        print("Motif\tFrequency")
        for g,h in v.items():
            print(g+"\t"+str(h))
        print("...\n")

elif q1=="7":
    z=uniq_g2
    print(v[:-2])
    motif_finder(mydict)
    for k,v in dict123.items():
        print(k)
        print("Motif\tFrequency")
        for g,h in v.items():
            print(g+"\t"+str(h))
        print("...\n")
else:
     print("error, please enter valid input")
     sys.exit()


