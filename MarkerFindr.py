import pandas as pd
import numpy as np
import argparse
import time

def get_argv() :
    parser = argparse.ArgumentParser(description='Compute the best 2 or 3 combination of genes or variants')
    parser.add_argument("-f", metavar='input', type=str, help='Variant or gene data in csv format')
    parser.add_argument("-m", metavar='metadata',type=str, help='Metadata in excel format')
    parser.add_argument("-o", metavar='output',type=str, help='Output name')
    parser.add_argument("-combination", metavar='metadata',type=str, help='Duo or trio combination')

    return parser.parse_args()

def make_addition(vec1,vec2) :
    return(np.add(vec1,vec2))

def make_soustraction(vec1,vec2) :
    res=[]
    for i in range(0,len(vec1)) :
        if vec1[i]==False and vec2[i]==False :
            res.append(False)
        elif vec1[i]==True and vec2[i]==True :
            res.append(False)
        else :
            res.append(True)
    return(res)


def compile_metadata(input, metadata) :
    data = pd.read_csv(input,sep ="\t")
    df = pd.read_excel(metadata)
    #df.columns = ["ID", "Cluster_name", "Source"]
    vect1=df.loc[df['Source'] == 'v1']
    vect2=df.loc[df['Source'] == 'v2']
    coln_vect1= vect1['ID'].unique()
    coln_vect2= vect2['ID'].unique()
    rowa=data['variant']
    res_vect1=pd.DataFrame(False,columns=coln_vect1, index=rowa)
    res_vect2=pd.DataFrame(False,columns=coln_vect2, index=rowa)
    for i in range(0,len(rowa)) :
            ids=data['k-samples'][i]
            ids_list=ids.split(",")
            for iden in ids_list :
                if iden in coln_vect1 :
                    res_vect1[iden][i]=True
                else :
                    res_vect2[iden][i]=True
    return res_vect1,res_vect2,rowa,coln_vect1,coln_vect2


def duo(res_vect1,res_vect2,rowa,output,coln_vect1,coln_vect2) :
    start_time = time.time()
    max_diff=0
    max_vec_vect1=""
    max_vec_vect2=""
    max_accuracy=0
    with open(output,"w") as fillout :
        for i in range(0,len(rowa)) :
            vec_vect1=res_vect1.values[i]
            vec_vect2=res_vect2.values[i]
            ##on skip si les gènes i sont présents partout
            if np.count_nonzero(vec_vect1)!=len(coln_vect1) and np.count_nonzero(vec_vect2)!=len(vec_vect2) :
                print("Gene numero " + str(i) + " sur " + str(len(rowa)))
                gene_i=rowa[i]
                for j in range(i, len(rowa)) :
                    vec_vect1=res_vect1.values[j]
                    vec_vect2=res_vect2.values[j]
                    ##on skip si les gènes j sont présents partout
                    if np.count_nonzero(vec_vect1)!=len(coln_vect1) and np.count_nonzero(vec_vect2)!=len(vec_vect2) :
                        gene_j=rowa[j]
                        if gene_i!=gene_j :
                            vec_vect1=make_addition(res_vect1.values[i], res_vect1.values[i])
                            if np.count_nonzero(vec_vect1)==0 :
                                vec_vect2=make_soustraction(res_vect2.values[i], res_vect2.values[i])
                            else :
                                vec_vect2=make_addition(res_vect2.values[i], res_vect2.values[i])
                        else :
                            vec_vect1=res_vect1.values[i]
                            vec_vect2=res_vect2.values[i]
                        nb_vect1=np.count_nonzero(vec_vect1)
                        nb_vect2=np.count_nonzero(vec_vect2)
                        accuracy=(nb_vect1+np.count_nonzero(vec_vect2==0) )/(len(vec_vect1)+len(vec_vect2))
                        line=gene_i+ "," + gene_j + "," + str(accuracy) + str(nb_vect1-nb_vect2)
                        if accuracy>max_accuracy :
                            max_diff=nb_vect1-nb_vect2
                            max_gene=gene_i+gene_j
                            max_vec_vect1=vec_vect1
                            max_vec_vect2=vec_vect2
                            max_accuracy=accuracy
                            line=gene_i+ ";" + gene_j + ";" + gene_k + ";" + str(accuracy)+ ";" + str(nb_vect1-nb_vect2) +";" + str(max_vec_vect1) + ";" + str(max_vec_vect2)
                            fillout.write(line +"\n")
    print("--- %s seconds ---" % (time.time() - start_time))

def trio(res_vect1,res_vect2,rowa,output,coln_vect1,coln_vect2) :
    start_time = time.time()
    max_diff=0
    max_vec_vect1=""
    max_vec_vect2=""
    max_accuracy=0
    with open(output,"w") as fillout :
        for i in range(0,len(rowa)) :
        #print("Gene numero " + str(i) + " sur " + str(len(rowa)))
            gene_i=rowa[i]
            vec_vect1=res_vect1.values[i]
            vec_vect2=res_vect2.values[i]
            if np.count_nonzero(vec_vect1)!=len(coln_vect1) and np.count_nonzero(vec_vect2)!=len(vec_vect2) :
                for j in range(i, len(rowa)) :
                    gene_j=rowa[j]
                    vec_vect1=res_vect1.values[j]
                    vec_vect2=res_vect2.values[j]
                    ##on skip si les gènes j sont présents partout
                    if np.count_nonzero(vec_vect1)!=len(coln_vect1) and np.count_nonzero(vec_vect2)!=len(vec_vect2) :
                        for k in range(j, len(rowa)) :
                            gene_k=rowa[k]
                            vec_vect1=res_vect1.values[k]
                            vec_vect2=res_vect2.values[k]
                            if np.count_nonzero(vec_vect1)!=len(coln_vect1) and np.count_nonzero(vec_vect2)!=len(vec_vect2) :
                                if gene_i==gene_j and gene_i!=gene_j :
                                    vec_vect1=make_addition(res_vect1.values[i], res_vect1.values[k])
                                    if np.count_nonzero(vec_vect1)==0 :
                                        vec_vect2=make_soustraction(res_vect2.values[i], res_vect2.values[k])
                                    else :
                                        vec_vect2=make_addition(res_vect2.values[i], res_vect2.values[k])

                                elif gene_i!=gene_j and gene_j==gene_k :
                                    vec_vect1=make_addition(res_vect1.values[i], res_vect1.values[k])
                                    if np.count_nonzero(vec_vect1)==0 :
                                        vec_vect2=make_soustraction(res_vect2.values[i], res_vect2.values[k])
                                    else :
                                        vec_vect2=make_addition(res_vect2.values[i], res_vect2.values[k])

                                elif gene_i==gene_k and gene_j!=gene_k :
                                    vec_vect1=make_addition(res_vect1.values[j], res_vect1.values[k])
                                    if np.count_nonzero(vec_vect1)==0 :
                                        vec_vect2=make_soustraction(res_vect2.values[j], res_vect2.values[k])
                                    else :
                                        vec_vect2=make_addition(res_vect2.values[j], res_vect2.values[k])

                                elif gene_i==gene_j and gene_i==gene_k :
                                    vec_vect1=res_vect1.values[i]
                                    vec_vect2=res_vect2.values[i]

                                else :
                                    vec_vect1=make_addition(res_vect1.values[i], res_vect1.values[j])
                                    vec_vect1=make_addition(vec_vect1, res_vect1.values[k])
                                    if np.count_nonzero(vec_vect1)==0 :
                                        vec_vect2=make_soustraction(res_vect2.values[i], res_vect2.values[j])
                                        vec_vect2=make_soustraction(vec_vect2, res_vect2.values[k])
                                    else :
                                        vec_vect2=make_addition(res_vect2.values[i], res_vect2.values[j])
                                        vec_pvect2=make_addition(vec_vect2, res_vect2.values[k])

                                nb_vect1=np.count_nonzero(vec_vect1)
                                nb_vect2=np.count_nonzero(vec_vect2)
                                accuracy=(nb_vect1+np.count_nonzero(vec_vect2==0) )/(len(vec_vect1)+len(vec_vect2))

                                if accuracy>max_accuracy :
                                    max_diff=nb_vect1-nb_vect2
                                    max_gene=gene_i+gene_j+gene_k
                                    max_vec_vect1=vec_vect1
                                    max_vec_vect2=vec_vect2
                                    line=gene_i+ ";" + gene_j + ";" + gene_k + ";" + str(accuracy)+ ";" + str(nb_vect1-nb_vect2) +";" + str(max_vec_vect1) + ";" + str(max_vec_vect2)
                                    max_accuracy=accuracy
                                    fillout.write(line +"\n")
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == "__main__":
    args = get_argv()
    res_vect1,res_vect2,rowa,coln_vect1,coln_vect2=compile_metadata(args.f, args.m)
    if args.combination=="duo" :
        duo(res_vect1,res_vect2,rowa,args.o,coln_vect1,coln_vect2)
    else :
        trio(res_vect1,res_vect2,rowa,args.o,coln_vect1,coln_vect2)
