import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage as ndimage
import sys
import cv2
import re
import Image

brightness = sys.argv[1]


alpha = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n']


def Interpolate(file1,file2):
   a = Image.open(file1)
   b = Image.open(file2)
   a = a.convert("RGBA")
   b = b.convert("RGBA")
   for i in range(1,10):
      new = Image.blend(a,b,i/10.0)
      print "Writing {0}{1}.png".format(file1.split(".")[0],alpha[i])
      new.save("{0}{1}.png".format(file1.split(".")[0],alpha[i]))



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
         return "output_{0}.png".format(number.zfill(8))

previous_file=""
thisFile = ""

for dat in sorted(os.listdir("output"), key= lambda x: int(("0" + x.split('.')[0].zfill(8)))) :
   print dat
   if dat.endswith(".dat"):
      try:

          this_file = convert_to_image(dat)
          if previous_file<>"":
             Interpolate(previous_file,this_file)

          previous_file = this_file

      except:
         pass



