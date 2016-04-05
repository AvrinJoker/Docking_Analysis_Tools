def overall_formating(l_rec,method,outfile,dire):
    '''format across different receptors to find the best methods in each case'''
    f = open(outfile,'w')
    f.write('Measure,%s\n'%method)
    for i in l_rec:
        l = [x for x in open('%s%s/%s_%s.csv'%(dire,i,i,method),'r')][2:]
        l = [x.split(',') for x in l]
        p = l[-1][1:]
        values = [(float(x[1]),x[2][:-1],x[0],l.index(x)) for x in l[:-1]]
        values = sorted(values,reverse=True)
        out = [values.index(x) for x in values if int(p[1]) in x or int(p[2]) in x]
        out =sorted(out)
        values = values[:out[1]]
        values = [','.join([x[2],str(x[0]),x[1]]) for x in values]
        values = ','.join(values)
        f.write('%s,%s\n'%(i,values))
    f.close()

if __name__=='__main__':
    overall_formating(['MAP4K4','ABL1'],'RMSD','test.csv','./rst/')