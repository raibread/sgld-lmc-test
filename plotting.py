from matplotlib import pyplot as plt
import numpy as np

L = 1e+6

# get results (seed 17 = dataset 1, seed 19 = dataset 2)

spath = "output/"

results = ( [spath + "results_seed17_reg1_n100.txt",
             spath + "results_seed17_reg2_n100.txt",
             spath + "results_seed17_reg3_n100.txt",
             spath + "results_seed19_reg1_n100.txt",
             spath + "results_seed19_reg2_n100.txt",
             spath + "results_seed19_reg3_n100.txt"] )
# get posterior mesh used in trapezoidal integration scheme             
grid = [spath + "PostGrid_seed17.txt", spath + "PostGrid_seed19.txt"]



# Plotting error bars for various samplers (SGLD, n = 1,2,4,8,16 and MHLD)
for k in [0,1,2,3,4,5]:
    fig, axs = plt.subplots(3, 2)
    for i in [0,1]:
        for j in [0,1,2]:
            # get results
            data = np.genfromtxt(results[i*3+j], delimiter = ",")
            
            ax = axs[j,i]
            if i == 0:
                ax.set_ylabel('R_'+str(j+1))
            if j == 0:
                ax.set_title('D_'+str(i+1))
            # indices of n = k results    
            ind = range(k*4+1,k*4+5)
            # estimate
            est = data[ind,4]
            # standard errors * 2
            err = 2*(data[ind,4] * (1 - data [ind,4]) / L * data[ind,-1])**.5
            # plot estimate and error
            ax.errorbar(range(1,5),est,yerr=err, fmt='s', ecolor='r', capthick = 2)
            # plot true statistic
            ax.plot([0,5], [data[0,0], data[0,0]],'k--')
            # plot reweighted statistic
            ax.plot(range(1,5),data[ind,5],'g+')
            
            if j == 2:
                if k != 5:
                    eps = data[ind,2:4]
                    labs = [ "["+str(x[0])+","+str(x[1])+"]" for x in eps ]
                else:
                    eps = data[ind,2]
                    labs = [ str(x) for x in eps ]
                ax.set_xticklabels(['']+labs+[''],rotation=-30,fontsize = 8)
            else:
                ax.set_xticklabels(['' for s in range(0,6)])

    plt.savefig('errors_n'+str(int(data[ind[0],1]))+'.pdf') 

# Test regions
R1 = np.array([[-.5, .25, .25, -.5, -.5],[.8,.8,1.8,1.8,.8]])
R2 = np.array([[ 0.  ,  1.  ,  1.  ,  0.  ,  0.  ],[-0.75, -0.75,  0.65,  0.65, -0.75]])
R3 = np.array([[ 0.9 ,  1.65,  1.65,  0.9 ,  0.9 ],[-1.9 , -1.9 , -0.9 , -0.9 , -1.9 ]])

xlimits = [(-1,2),(-1,2)]
ylimits = [(-2,2),(-2.5,2.5)]


# Plotting posterior contour plots for both datasets
# region boxes overlayed
fig, axs = plt.subplots(1,2)
for l in [0,1]:
    ax = axs[l]
    # load posterior mesh
    data = np.genfromtxt(grid[l], delimiter = ",")
    # plot posterior mesh used to integrate true statistic
    ax.contour(data[:,0].reshape([999,999])[:,0],
                data[:,1].reshape([999,999])[0,:],
                data[:,2].reshape([999,999]).T) 
    # define regions over which we estimate probability
    r1, = ax.plot(R1[0,:],R1[1,:])
    r2, = ax.plot(R2[0,:],R2[1,:]) 
    r3, = ax.plot(R3[0,:],R3[1,:])
    ax.legend([r1,r2,r3],['R1','R2','R3'])
    ax.set_title('Dataset '+str(l+1))

    ax.set_xlim(xlimits[l])
    ax.set_ylim(ylimits[l])
plt.savefig('distr_dataset.pdf')

