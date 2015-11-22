import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage as ndimage
import sys
import cv2








for file in os.listdir("output"):
   if file.endswith(".dat"):
      dat = file









      path = "output/" + file
      print path

      x = np.loadtxt(path ,skiprows=1)
      with open(path,'r') as f:
           w,h = f.readline().split()



      w,h = int(w),int(h)
      plt.figure(figsize=(w/20,h/20))

      y = np.reshape(x, (h,w,3)) * 255


      plt.imshow(y,origin='lower',extent=(0,w,0,h))
      plt.savefig("output_{0}.png".format(file))
      plt.close()



