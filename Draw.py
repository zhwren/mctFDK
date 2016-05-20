############################################################
# Author       : ZhuHaiWen                                 #
# Email        : zhuhw@ihep.ac.cn/zhwren0211@whu.edu.cn    #
# Last modified: 2015-12-03 8:48:1449132507
# Filename     : Draw.py
# Phone Number : 18625272373                               #
# Discription  :                                           #
############################################################

import struct
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

fmt = str(512*512*10) + 'f'
class mct:
    """
    Plot the reconstructed data.
    """
    def __init__(self,plotNumber=8,fileName='ReconData.rcn'):
        self.plotNumber = plotNumber
        self.imgData    = np.zeros([512,512])
        self.dataFile   = open(fileName,'rb')
        self.data       = struct.unpack(fmt,self.dataFile.read(4*512*512*10))
        self.dataFile.close()

    def read(self,n=-1):
        for i in range(512):
            for j in range(512):
                if n == -1:
                    n = self.plotNumber -1
                self.imgData[i,j] = self.data[i*512+j+512*512*n]

    def draw(self):
        if self.plotNumber != 0:
            self.read()
            cmap = mpl.cm.gray
            norm = mpl.colors.Normalize(vmin=0.019,vmax=0.024)
            plt.figure(figsize=(6,6))
            plt.imshow(self.imgData,cmap=cmap,norm=norm)
#            plt.figure(figsize=(6,6))
#            plt.plot(range(512), self.imgData[256,:])
            plt.show()
        else:
            for i in range(10):
                print "picture number:",(i+1)
                self.read(i)
                cmap = mpl.cm.gray
                norm = mpl.colors.Normalize(vmin=0.019,vmax=0.023)
                plt.imshow(self.imgData,cmap=cmap,norm=norm)
                plt.show()

    def add(self):
        for n in range(10):
            for i in range(512):
                for j in range(512):
                    self.imgData[i,j] = self.imgData[i,j]+self.data[i*512+j+512*512*n]/10
        cmap = mpl.cm.gray
        norm = mpl.colors.Normalize(vmin=0.019,vmax=0.023)
        plt.figure(figsize=(6,6))
        plt.imshow(self.imgData,cmap=cmap,norm=norm)
#        plt.figure(figsize=(6,6))
#        plt.plot(range(512), self.imgData[256,:])
        plt.show()

if __name__ == "__main__":
    mct(0).draw()
