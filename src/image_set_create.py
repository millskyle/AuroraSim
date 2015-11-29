import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage as ndimage
import sys
import cv2
import re

brightness = sys.argv[1]






def convert_to_image(dat):
         print dat
         path = "output/" + dat
         x = np.loadtxt(path ,skiprows=1)
         number = re.findall("\d*",dat)[0]
         with open(path,'r') as f:
            w,h = f.readline().split()

         w,h = int(w),int(h)
         plt.figure(figsize=(w/20,h/20))

         y = np.reshape(x, (h,w,3)) * 255
         y = y*float(brightness)

         plt.imshow(y,origin='lower',extent=(0,w,0,h))
         plt.savefig("output_{0}.png".format(number.zfill(8)))
         plt.close()



for dat in os.listdir("output"):
   if dat.endswith(".dat"):
      try:
         convert_to_image(dat)
      except:
         pass



