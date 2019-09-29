def get_expression(file_name):
    fh = open(file_name, 'r')
    # You will get the column head labels and save them into the following variables
    col1, col2, col3 = fh.readline().strip().split('\t')
    # You need one more column label
    col4 = 'Exp'
    dict0 = {}
    lines = fh.readlines()
    for i in lines:
        x=i.strip().split('\t')
        # You can determine the expression status such as +, -, and other cases such as "s"
        if int(x[1])/int(x[2]) >= 1.5:
            exp = '+'
        elif int(x[2])/int(x[1]) >= 1.5:
            exp = '-'
        else:
            exp = '.'
        # You insert the concatenated string as the dictionary value
        dict0[x[0]]=x[1]+':'+x[2]+':'+exp
    return dict0, col1, col2, col3, col4
get_expression('cancer_lung_expression.txt')