data = np.array([images[i,:,:].flatten() for i in range(images.shape[0])])
avgim = np.mean(data,axis=0)
residual = np.array([data[i,:]-avgim for i in range(data.shape[0])])
COVMAT = np.cov(residual.T)
vals,vecs = np.linalg.eig(COVMAT) # this peice takes a long time.

Then remove mean, then remove a handfull of np.dot(eigvecs,data) and plot the >0 pixels
