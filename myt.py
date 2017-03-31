import string
import datetime

#count Protein_coding genes from raw-reading gff file. just checks for presence of the gene field and gene_biotype=Protein_coding
#does not parse the file for a faster count
#uses a file handle as input
def count_protein_coding(handle):
    counter = 0
    for line in handle.readlines():
        if ("\tgene\t" in line) and ("gene_biotype=Protein_coding" in line):
            counter +=1
    handle.seek(0)
    return (counter)

#%%

#generates a new GFF file with only gene field and gene_biotype=Protein_coding positive strings
#does not parse the file
def extract_protein_coding(input_handle):
    pc = open ("protein_coding.gff3", "w")
    pc.write("##gff-version 3\n")
    for line in input_handle.readlines():
        if ("\tgene\t" in line) and ("gene_biotype=Protein_coding" in line):
            pc.write(line)
    pc.close()
    input_handle.seek(0)
#%%

#parses the file created by extract_protein_coding()
def parse_p_c():
    dic = {}
    v = [] #line-parsing vector
    handle = open("protein_coding.gff3", "r")
    print(handle.readline())
    for line in handle.readlines():
        v = str(line).split('\t')
        print(v[0], end='\r')
        scaffold = str(v[0])
        gene = str(v[8]).split(";")
        gene = str(gene[0]).strip("ID=")
        if (scaffold not in dic.keys()):
            #dic[v[0]] = [str(str(v[8]).split(";")[0]).strip("ID=")]
            dic[scaffold] = []
            tmp = list(dic[scaffold])
            tmp.append(gene)
            dic[scaffold] = tmp
        else:
            tmp = list(dic[scaffold])
            tmp.append(gene)
            dic[scaffold] = tmp
            #dic.setdefault(v[0],[]).append(str(str(str(v[8]).split(";")[0]).strip("ID=")))

    c = 0

    for i in dic:
        print("%s : %s" %(i, dic[i]))
        c = c+1

    print("\nscaffolds totali : %d" %c)
    handle.close()
    return dic

#%%

if __name__ == '__main__':
    in_scaffold = open("MGAL10A.annot.gff3", "r")
    print(count_protein_coding(in_scaffold))
    extract_protein_coding(in_scaffold)
    s = parse_p_c()
    d = {}
    PAV = open("PAV_list.txt", "r")
    raw_output = open("raw_output.csv","w")
    output = open("output.csv", "w")
    for gene in PAV.readlines():
        gene = gene.strip()
        found = False
        for i in s:
            if gene in s[i]:
                found = True
                scaffold = i
                raw_output.write("%s,%s\n" %(gene, scaffold))
                if scaffold not in d.keys():
                    d[scaffold] = []
                    tmp = list(d[scaffold])
                    tmp.append(gene)
                    d[scaffold] = tmp
                else:
                    tmp = list(d[scaffold])
                    tmp.append(gene)
                    d[scaffold] = tmp
                #print("trovato %s in %s" %(gene, i))
        if not found:
            print ("\nnon ho trovato %s\n", gene)
        else:
            print ("trovato %s" %gene,end='\r')

        #print("%s :: %s" %(i, s[i]))

    #print (list(s.values()))
    print("\ncreazione file")
    output.write("scaffold,geni_totali,geni_trovati\n")
    for scaffold in s:
        totali=len(s[scaffold])
        if scaffold in d:
            trovati=len(d[scaffold])
            print("%s\t%d\t%d" %(scaffold,totali,trovati))
            output.write("%s,%d,%d\n" %(scaffold,totali,trovati))
        else:
            trovati=0
            print("%s\t%d\t%d" %(scaffold,totali,trovati))
            output.write("%s,%d,%d\n" %(scaffold,totali,trovati))



    PAV.close()
    in_scaffold.close()
    raw_output.close()
    output.close()
#%%
