vls3=vlsdata[:,:1<<8]
vls2=vlsdata[:,((1<<10)+(1<<9)):((1<<11)-(1<<8))]
covmat=np.cov(vls2.T,vls3.T)
ridge = [np.argmax(covmat[i,:]) for i in range(covmat.shape[0])]
plt.plot(ridge,'.');plt.show()
