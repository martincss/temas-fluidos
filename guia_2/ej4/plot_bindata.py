import numpy as np
import matplotlib.pyplot as plt

# Reads a binary file and plots a cut in the x-y plane.
# Execute in ipython with '%run plot_bindata.py'

# Path to the binary data
path = './N_10/'

# Spatial resolution
NX = 128
NY = 128
NZ = 64
shape = (NX,NY,NZ)

# Reads binary files
#vx = np.fromfile(path+'vx.0001.out',dtype=np.float32).reshape(shape,order='F')

# Show a horizontal cut of the field in the middle of the box
# plt.figure(1)
# plt.imshow(vx[:,:,NZ/2])
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()
plt.figure()
for i in range(1,9):
    vx = np.fromfile(path+'vz.{:04d}.out'.format(i),dtype=np.float32).reshape(shape,order='F')
    plt.clf()
    plt.title('i = {:}'.format(i))
    plt.imshow(vx[:,NY//2,:], cmap = 'jet')#, vmin = , vmax = )
    plt.xlabel('x')
    plt.ylabel('y')
    plt.colorbar()
    plt.pause(1)
    print(np.max(vx))
plt.show()
