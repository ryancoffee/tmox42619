>>> fig,axs = plt.subplots(7,1,sharex=True,sharey=True)
>>> for g in range(1,8):
...     axs[g-1].imshow(np.array(data[g])[:,100:],vmin=0,vmax=5e-8,origin='lower')
... 
<matplotlib.image.AxesImage object at 0x7f613c703c10>
<matplotlib.image.AxesImage object at 0x7f613c77ab50>
<matplotlib.image.AxesImage object at 0x7f613c77a850>
<matplotlib.image.AxesImage object at 0x7f613c755370>
<matplotlib.image.AxesImage object at 0x7f613c77a2b0>
<matplotlib.image.AxesImage object at 0x7f613c75a070>
<matplotlib.image.AxesImage object at 0x7f613c75a370>
>>> plt.show()
 data = [ [f['port_12']['hist'][()][i,g,:]/(f['vls']['norm'][()][i]*(f['port_12']['qbins'][()][1:] - f['port_12']['qbins'][()][:-1])) for i in range(f['vls']['norm'].shape[0])] for g in range(8)]
