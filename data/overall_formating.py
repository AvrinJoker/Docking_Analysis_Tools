import os
def overall_formating(l_rec,method,outfile,best=True, meth=None, dire='.'):
    '''format across different receptors to find the best methods in each case or summarize a specific method'''
    f = open(outfile,'w')
    f.write('Measure,%s\n'%method)
    for i in l_rec:
        l = [x for x in open('%s%s_%s.csv'%(dire,i,method),'r')][2:]
        l = [x.split(',') for x in l]
        p = l[-1][1:]
        if best == True:
            values = [(float(x[1]),x[2][:-1],x[0],l.index(x)) for x in l[:-1]]
            values = sorted(values,reverse=True)
            out = [values.index(x) for x in values if int(p[1]) in x or int(p[2]) in x]
            out =sorted(out)
            values = values[:out[1]]
            values = [','.join([x[2],str(x[0]),x[1]]) for x in values]
            values = ','.join(values)
            f.write('%s,%s\n'%(i,values))
        else:
            values = [','.join(x) for x in l[:-1] if meth in x][0]
            print values
            f.write('%s,%s'%(i,values))
    f.close()

if __name__=='__main__':
    l_names = [x.split('_')[0] for x in os.listdir('.') if x.endswith('csv')]
    overall_formating(l_names,'IC50_R2','./figure/overall.csv')
    for i in ['min-cross','align-cross','dock-cross','align-close','dock-close']:
        overall_formating(l_names,'IC50_R2','./figure/%s.csv'%i,False, i)
