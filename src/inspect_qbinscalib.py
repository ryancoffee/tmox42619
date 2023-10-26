>>> for p in portkeys:
...     qmeans[p] = (q[p]['qbins'][()][:-1]+q[p]['qbins'][()][1:])*0.5
...     qdiffs[p] = q[p]['qbins'][()][1:]-q[p]['qbins'][()][:-1]
...     plt.plot(qmeans[p],1./qdiffs[p],label=p)
... 
>>> plt.legend()
<matplotlib.legend.Legend object at 0x7f1fe54b5d30>
>>> plt.show()
>>> portkeys = [k for k in q.keys() if re.search('port',k)]
qname = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output_multiretardation/h5files/quantHist_nonuniform_runs132-133-134.NNO.qknob0.000.tmox42619.h5'
>>> q = h5py.File(qname,'r')
>>> q.keys()

