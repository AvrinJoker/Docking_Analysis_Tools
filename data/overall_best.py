import matplotlib.pyplot as plt
import numpy as np

def infile_parser(infile):
    lib = {'RMSD':'Percentage of RMSD <= 2.0 A','Kd_SPM':'Kd Spearman Correlation', 'IC50_SPM':'IC50 Spearman Correlation', 'Kd_R2':'Kd R2', 'IC50_R2':'IC50 R2'}
    l = [x.split(',') for x in open(infile,'r')]
    measure_name = lib[l[0][1][:-1]]
    
    value = [(float(x[2]),x[0],float(x[3])) for x in l[1:]]
    value = sorted(value, reverse=True)
    rec_name = [x[1] for x in value]
    y_value = [float(x[0]) for x in value]
    y_value = [float('%.2f'%x) for x in y_value]
    y_value = np.asarray(y_value)
    y_error = [float(x[2]) for x in value]
    y_error = np.asarray(y_error)
    return measure_name,rec_name,y_value,y_error

def overall_best(infile, outfile):
    '''plotting best scores across different receptors from high to low

    data includes RMSD, Spearman (IC50 & Kd) and R2 (IC50 & Kd)
    '''
    
    # get plotting data from file
    measure_name,rec_name,y_value,y_error = infile_parser(infile)
    figsize = ((len(y_value)+1)*0.5,3)

    # set figure basic info
    fig = plt.figure(figsize=figsize, dpi=600)
    ax = fig.add_subplot(111,ylim=(0,1),xlim=(-1,len(y_value)))
    x = [x for x in range(len(y_value))]

    # add error bar
    ax.errorbar(x,y_value, yerr=y_error, ecolor= 'k', capthick=2, fmt=None, capsize=3, elinewidth=2, markeredgewidth=1)
    ax.scatter(x, y_value, facecolor='k', s=100, linewidths=0, alpha=0.8)
    
    # set the axises & format the fig
    ax.set_ylabel(measure_name,fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(rec_name)
    plt.xticks(rotation='vertical')
    
    # mark the x value
    for i in range(len(y_value)):
        ax.annotate(y_value[i], xy=(i,y_value[i]-y_error[i]-0.07),ha='center',fontsize=8)

    # draw a 0.5 line
    ax.axvspan(-1, len(y_value), ymin=0.5, ymax=0.5, alpha=0.3, linestyle='dotted', linewidth=2)


    # save figure
    plt.savefig(outfile,dpi=600,bbox_inches='tight')

def all_in_one(outfile,dire):
    # set figure basic info
    measure_name,rec_name,y_value,y_error = infile_parser('%soverall.csv'%dire)
    
    figsize = ((len(y_value)+1)*0.15,(len(y_value)+1)*0.15)
    fig = plt.figure(figsize=figsize, dpi=600)
    ax = fig.add_subplot(111,ylim=(0,1),xlim=(0,1))
    x = [x/float(len(y_value)-1) for x in range(len(y_value))]
    ax.set_ylabel(measure_name,fontweight='bold', fontsize=16)
    ax.set_xlabel('Protein percentage', fontweight='bold', fontsize=16)
    ax.axvspan(0, 1, ymin=0.5, ymax=0.5, alpha=0.3, linewidth=2)
    cspace = ['#6495ed','#1e90ff','#4169e1','#c71585','#db7093']
    nspace = ['min-cross','align-cross','dock-cross','align-close','dock-close']
    #overall
    ax.fill_between(x,y_value-y_error,y_value+y_error, facecolor = 'k', linewidths = 0, alpha = 0.2)
    ax.plot(x, y_value, lw=2, alpha=0.8, color='k', label = 'overall')
    ax.scatter(x, y_value, facecolor='k', s=10, linewidths=0, alpha=0.8)
    
    w = [i for i in y_value if i >= 0.5]
    w = float(len(w)-0.5)/(len(y_value)-1)
    ax.axvline(w, ymin=0, ymax=0.5, alpha=0.3, color= 'k', linewidth=1 )
    ax.annotate('%.2f'%w, xy=(w,0.02),ha='center',fontsize=8)

    
    for i in range(len(nspace)):
        cc = cspace[i]
        measure_name,rec_name,y_value,y_error = infile_parser('%s%s.csv'%(dire,nspace[i]))
        ax.fill_between(x,y_value-y_error,y_value+y_error, facecolor = cc, linewidths = 0, alpha = 0.2)
        ax.plot(x, y_value, lw=2, alpha=0.8, color=cc, label='%s'%nspace[i])
        ax.scatter(x, y_value, facecolor=cc, s=10, linewidths=0, alpha=0.8)
        w = [j for j in y_value if j >= 0.5]
        w = float(len(w)-0.5)/(len(y_value)-1)
        ax.axvline(w, ymin=0, ymax=0.5, alpha=0.3, color= 'k', linewidth=1, )
        ax.annotate('%.2f'%w, xy=(w,0.02+0.01*i),ha='center',fontsize=8)

    plt.legend(fontsize=8)
    plt.savefig(outfile,dpi=600,bbox_inches='tight')
    plt.close()


if __name__=='__main__':
    #overall_best('./figure/overall.csv', 'overall.png')
    #for i in ['min-cross','align-cross','dock-cross','align-close','dock-close']:
    #    overall_best('./figure/%s.csv'%i, '%s.png'%i)
    all_in_one('IC50_R2_methods_all_compared.png','./figure/')
