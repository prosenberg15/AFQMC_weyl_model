import numpy as np
import matplotlib.pyplot as plt
import sys

Ns=8

data=np.array(np.loadtxt(sys.argv[1]))
#data2=np.array(np.loadtxt(sys.argv[2]))

asites=np.arange(1,2*Ns,2)
bsites=np.arange(2,2*Ns+1,2)

#ni^A nj^A
plt.errorbar(asites,2.0*(data[:Ns,1]+data[4*Ns:Ns+4*Ns,1]),yerr=data[:Ns,3])
#ni^A nj^B
plt.errorbar(bsites,2.0*(data[Ns:2*Ns,1]+data[Ns+4*Ns:2*Ns+4*Ns,1]),yerr=data[:Ns,3])


#plt.errorbar(data[:Ns,0],(data[:Ns,1]+data[4*Ns:Ns+4*Ns,1]+data[Ns:2*Ns,1]+data[Ns+4*Ns:2*Ns+4*Ns,1]+data[2*Ns:3*Ns,1]+data[2*Ns+4*Ns:3*Ns+4*Ns,1]+data[3*Ns:4*Ns,1]+data[3*Ns+4*Ns:4*Ns+4*Ns,1]),yerr=data[Ns+4*Ns:2*Ns+4*Ns,3])
#plt.errorbar(data2[:Ns,0],2.0*(data2[:Ns,1]+data2[4*Ns:Ns+4*Ns,1]+data2[Ns:2*Ns,1]+data2[Ns+4*Ns:2*Ns+4*Ns,1]+data2[2*Ns:3*Ns,1]+data2[2*Ns+4*Ns:3*Ns+4*Ns,1]+data2[3*Ns:4*Ns,1]+data2[3*Ns+4*Ns:4*Ns+4*Ns,1]),yerr=data2[Ns+4*Ns:2*Ns+4*Ns,3])

#plt.plot(np.arange(0,Ns+1,0.001),len(np.arange(0,Ns+1,0.001))*[1.0])
#plt.xlim(1,Ns)
#plt.ylim(1.9,2.1)
#plt.xticks(np.arange(1,33,1))
plt.xticks(np.arange(1,2*Ns+1,1))
#plt.ylim(0.2,0.3)

plt.show()
