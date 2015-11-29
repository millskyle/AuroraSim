import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage as ndimage
import sys
import cv2


dat = sys.argv[1]
alter = sys.argv[2]

path = dat
print path

x = np.loadtxt(path ,skiprows=1)
with open(path,'r') as f:
     w,h = f.readline().split()



w,h = int(w),int(h)
plt.figure(figsize=(w/20,h/20))

y = np.reshape(x, (h,w,3)) * 255
y = y*float(alter)

y[np.where(y > 255)] = 0

plt.imshow(y,origin='lower',extent=(0,w,0,h))
plt.savefig("output.png")
plt.close()



