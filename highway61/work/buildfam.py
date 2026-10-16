import sys,pickle,mkfont
WX={'B':4,'C':4,'D':4,'E':5,'F':5,'EM':0,'ER':0,'EMC':0}
SB={'B':(0.72,0.83),'C':(0.85,0.94),'D':(1.0,1.03),'E':(1.04,1.05),'F':(1.06,1.08),'EM':(1.04,1.05),'ER':(1.04,1.05),'EMC':(1.04,1.05)}
s=sys.argv[1]
import json,string
DL=json.load(open('spacing_delta.json')).get(s)
EX={}
if DL:
    for c in string.ascii_lowercase: EX[c]=DL['lc']
    for c in string.ascii_uppercase: EX[c]=DL['uc']
    for c in '0123456789': EX[c]=DL['fig']
EX['w']=EX.get('w',0)+WX[s]
if s=='B': EX['x']=EX.get('x',0)+4
mkfont.EXTRA_SB=EX
mkfont.SMW=5 if s in ('B','C','D','F') else 7
src=dict(pickle.load(open('caps_%s.pkl'%({'ER':'E','EMC':'EM'}.get(s,s)),'rb'))); src.update(pickle.load(open('lcE_%s.pkl'%s,'rb')))
import numpy as np,cv2
bi,ox,oy,t=src['i']; n,lab,st,_=cv2.connectedComponentsWithStats(bi)
k=min(range(1,n),key=lambda j:st[j][1]); dot=(lab==k).astype(np.uint8)
src['.']=(dot,ox,float(st[k][1]+st[k][3]),'period')
# apostrophe: top of the Y's right arm, cut flat (horizontally) top and bottom
yb,yox,yoy,_=src['V']
r0=int(round(yoy-690)); r1=int(round(yoy-470))
crop=np.zeros_like(yb); crop[max(0,r0-2):r1+1]=yb[max(0,r0-2):r1+1]
n,lab,st,cen=cv2.connectedComponentsWithStats(crop)
k=max(range(1,n),key=lambda j:cen[j][0]); arm=(lab==k).astype(np.uint8)
hb,hox,hoy,_=src['H']
stem=float(np.median([len(r) for r in np.split(np.nonzero(hb[int(hoy-150)])[0],np.where(np.diff(np.nonzero(hb[int(hoy-150)])[0])!=1)[0]+1)]))
h=int(round(1.7*stem))
ap=np.zeros_like(yb); rb=int(round(yoy-(690-h)))
ap[r0-2:rb]=arm[r0-2:rb]
xs=np.nonzero(ap.any(0))[0]; ap=ap[:,xs[0]:xs[-1]+1]
midrow=np.nonzero(ap[int(round(yoy-(690-h/2)))])[0]; run=midrow[-1]-midrow[0]+1
ap=np.pad(ap,((0,0),(20,20)))
for ch in ["'",'\u2019','\u2018']: src[ch]=(ap,20,yoy,'apos')
shift=int(round(run+0.55*stem)); dq=np.zeros((ap.shape[0],ap.shape[1]+shift),np.uint8); dq[:,:ap.shape[1]]=ap; dq[:,shift:shift+ap.shape[1]]|=ap
for ch in ['"','\u201d','\u201c']: src[ch]=(dq,20,yoy,'quotes')
aph=yoy-np.nonzero(ap.any(1))[0].max()   # y of apostrophe bottom
src[',']=(ap,20,yoy-(690-125),'comma')
# hyphen: the H crossbar between the stems, centred on half the x-height
r=hb[int(hoy-345)]
if True:
    rows=[y for y in range(200,500) if (np.diff(np.nonzero(hb[int(round(hoy-y))])[0])!=1).sum()==0]
    yb0,yb1=min(rows),max(rows)
    sr=hb[int(round(hoy-150))]; xs2=np.nonzero(sr)[0]; sp=np.where(np.diff(xs2)!=1)[0]
    xl=xs2[sp[0]]+1; xr=xs2[sp[0]+1]-1
    hy=np.zeros((yb1-yb0+1,xr-xl+1),np.uint8)+1
    hy=np.pad(hy,((10,10),(10,10)))
    src['-']=(hy,10,10+(yb1-yb0+1)+262-(yb1-yb0+1)/2-0.0,'hyphen')
# exclamation: the i turned upside down
ytop=oy-np.nonzero(bi.any(1))[0].min(); H=bi.shape[0]
src['!']=(np.ascontiguousarray(bi[::-1,::-1]),ox,ytop-oy+H-1,'exclam')
cs,ls=SB[s]
OUT=('/mnt/user-data/outputs/Highway1966-EM-Constructed.ttf','Highway 1966 E(M) Constructed') if s=='EMC' else ('/mnt/user-data/outputs/Highway1966-E-Reference.ttf','Highway 1966 E Reference') if s=='ER' else ('/mnt/user-data/outputs/Highway1966-EM.ttf','Highway 1966 E(M)') if s=='EM' else ('/mnt/user-data/outputs/HighwayEven1966-%s.ttf'%s,'Highway Even 1966 %s'%s)
mkfont.build(OUT[0],family=OUT[1],src=src,sbscale=cs,lcscale=ls,
  lcadd=0,desc=('Series E(M) capitals and figures traced from the 1983 Standard Alphabets, with the original 1966 lower case at its drawn weight, and punctuation fabricated from the letters.' if s=='EM' else None) or '1966 FHWA Series %s capitals and figures traced from the 1966 Standard Alphabets, with an even-stroke lowercase derived from the 1966 lower case by centreline reconstruction.'%s)
print(s,'done')
