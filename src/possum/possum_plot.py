#!/usr/bin/env fslpython

import numpy as np
import sys
# also tries to import malplotlib and PIL below

# parse arguments
if len(sys.argv)<=2:
    print("Usage: %s <input file> <output basename>" % sys.argv[0])
    sys.exit(0)
 
filenm=sys.argv[1]
tmpnam=sys.argv[2]

mplot=True
try:
    import matplotlib.pyplot as plt
except ImportError:
    print('Failed matplotlib import: please install matplotlib to get improved plots')
    mplot=False
try:
    from PIL import Image
except ImportError:
    print('Failed PIL import: please install Pillow to get improved plots')
    mplot=False
    
if mplot:
    #print('Found matplotlib and PIL')
    m=np.loadtxt(filenm)
    fig=plt.figure()
    fig.set_dpi(100)
    fig.set_figheight(3)
    fig.set_figwidth(5)
    plt.plot(m[:,0],m[:,1:])
    plt.xlabel('Time',fontsize=20)
    plt.ylabel('Magnitude',fontsize=20)
    plt.savefig(tmpnam+".png",bbox_inches="tight",pad_inches=0.2)
    im=Image.open(tmpnam+".png")
    im.save(tmpnam+".gif")
else:
    #print('Cannot find matplotlib or PIL: reverting to fsl_tsplot')
    import os
    fsldir=os.environ.get('FSLDIR')
    os.system(fsldir+"/bin/fsl_tsplot -i "+filenm+" -o "+tmpnam+".png --start=2 --finish=7 -x MotionMatrixRow -y Magnitude -w 500")
    os.system(fsldir+"/bin/pngappend "+tmpnam+".png "+tmpnam+".gif")

