import sys,hroverlay,hires
s=sys.argv[1]
P={'B':['0004','0005','0006','0007'],'C':['0008','0009','0010','0011'],'D':['0012','0013','0014','0015'],'F':['0020','0021','0022','0023']}[s]
rows=[['ABC','DEF','GHI'],['JKL','MNO','PQR'],['STU','VWX','YZ1'],['234','567','890']]
for p,r in zip(P,rows): hires.LAY[p]=r
hroverlay.LAY=hires.LAY
import pickle,numpy as np
hroverlay.font='/mnt/user-data/outputs/HighwayEven1966-%s.ttf'%s
hroverlay.overlay(P,690,float(sys.argv[2]),'/mnt/user-data/outputs/hires-capitals-%s-vs-scan.png'%s)
