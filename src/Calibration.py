#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.5.7-py39/bin/python3

import numpy as np
import matplotlib.pyplot as plt

class CalibData:
    def __init__(self):
        self.portkeys = ['port_0','port_1','port_12','port_13','port_14','port_15','port_4','port_5']
        self.data = {}
        self.data['port_0'] = {}
        self.data['port_0']['q'] = np.array([126,150,180,10.8,23.3,28.7])
        self.data['port_0']['e'] = np.array([371.5,368.5,365,502.5,496,490])
        self.data['port_1'] = {}
        self.data['port_1']['q'] = np.array([138.6,171.5,195,18.5,31.5,41.27,7.9,4.8])
        self.data['port_1']['e'] = np.array([371.5,368.5,365.5,502.5,496,490,558.4,562.6])
        self.data['port_12'] = {}
        self.data['port_12']['q'] = np.array([160,136,71.5,56.8,50.7,24,19.3,14.1,6.8,2.4])
        self.data['port_12']['e'] = np.array([187.45,191.2,365,368.5,371,490,496,502.5,558.4,562.6])
        self.data['port_13'] = {}
        self.data['port_13']['q'] = np.array([120,55.2,48.7,30.7,24.5,22,9.6,7,5.1,3,2])
        self.data['port_13']['e'] = np.array([59,187.45,191.2,365,368.5,371,490,496,502.5,558.4,562.6])
        self.data['port_14'] = {}
        self.data['port_14']['q'] = np.array([133,87.9,59.3,28.3,12])
        self.data['port_14']['e'] = np.array([490,496,502.5,558.4,562.6])
        self.data['port_15'] = {}
        self.data['port_15']['q'] = np.array([93,55.3,35.3,16.7, 10.8])
        self.data['port_15']['e'] = np.array([490,496,502.5,558.4,562.6])
        self.data['port_4'] = {}
        self.data['port_4']['q'] = np.array([220,202.7,112,91.4,80,40,29.1,20.6,9.3,2.7])
        self.data['port_4']['e'] = np.array([187.45,191.2,365,368.5,371,490,496,502.5,558.4,562.6])
        self.data['port_5'] = {}
        self.data['port_5']['q'] = np.array([164,68.6,60.1,38.5,29.8,26.7,12,9.1,6.3,3,1])
        self.data['port_5']['e'] = np.array([59,187.45,191.2,365,368.5,371,490,496,502.5,558.4,562.6])

        self.wins = {'o1s':[0,100],'n1s':[100,300],'nam':[300,400],'oam':[400,525],'IV':[525,600]}
        self.theta = {}
        for k in self.portkeys:
            self.theta[k] = {}

    def fit(self):
        for k in self.portkeys:
            for w in self.wins.keys():
                x = self.data[k]['q']
                y = self.data[k]['e'] 
                inds = np.where((y>self.wins[w][0])*(y<self.wins[w][1]))
                X = np.ones((2,len(inds[0])))
                Y = y[inds]
                if X.shape[0] > 1:
                    for i in range(X.shape[0]):
                        X[i,:] = np.power(x[inds],int(i))
                    self.theta[k][w] = np.dot(Y,np.linalg.pinv(X))
                    printstring = '%s fit\t%s (%i,%i)\t:\t'%(k,w,self.wins[w][0],self.wins[w][1])
                    for v in self.theta[k][w]:
                        printstring += '%.3f\t'%v
                    print(printstring)
        return self


def main():
    calib = CalibData()
    plotting = False
    if plotting:
        print('suspect O2s and N2s at 558.4 and 562.6 respectively')
        _=[ plt.plot(calib.data[k]['q'],calib.data[k]['e'],'o',label=k) for k in calib.portkeys[0:8:2] ]
        plt.legend()
        plt.show()
    calib.fit()
    return

if __name__ == '__main__':
    main()
