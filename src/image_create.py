import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage as ndimage
import sys
import cv2

dat = sys.argv[1]

x = np.loadtxt(dat ,skiprows=1)

with open(dat,'r') as f:
     w,h = f.readline().split()



w,h = int(w),int(h)
plt.figure(figsize=(w/20,h/20))

y = np.reshape(x, (h,w,3)) * 255



#y[:,:,0] = ndimage.gaussian_filter(y[:,:,0], 10)
#y[:,:,1] = ndimage.gaussian_filter(y[:,:,1], 10)
#y[:,:,2] = ndimage.gaussian_filter(y[:,:,2], 10)

#y[:,:,0] = cv2.blur(y[:,:,0], (11,11))
#y[:,:,0] = cv2.blur(y[:,:,0], (11,11))
y = cv2.blur(y, (2,2))


plt.imshow(y,origin='lower',extent=(0,w,0,h))
plt.savefig('output.png')



