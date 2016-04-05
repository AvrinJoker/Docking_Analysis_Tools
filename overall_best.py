import matplotlib.pyplot as plt

def overall_best(infile, outfile):
    '''plotting best scores across different receptors from high to low

    data includes RMSD, Spearman (IC50 & Kd) and R2 (IC50 & Kd)
    '''
    def _infile_parser(infile):
        lib = {'RMSD':'Percentage of RMSD <= 2.0 A','SPM_Kd':'Kd Spearman Correlation', 'SPM_IC50':'IC50 Spearman Correlation', 'R2_Kd':'Kd R2', 'R2_IC50':'IC50 R2'}
        l = [x.split(',') for x in open(infile,'r')]
        measure_name = lib[l[0][1][:-1]]
        
        value = [(x[2],x[0],x[3]) for x in l[1:]]
        value = sorted(value, reverse=True)
        rec_name = [x[1] for x in value]
        y_value = [float(x[0]) for x in value]
        y_error = [float(x[2]) for x in value]
        return measure_name,rec_name,y_value,y_error
    
    # get plotting data from file
    measure_name,rec_name,y_value,y_error = _infile_parser(infile)
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


if __name__=='__main__':
    overall_best('test.csv', 'overall.png')
