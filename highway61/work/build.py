import pickle,numpy as np,cv2
from src import *
L=pickle.load(open('L.pkl','rb'))
G=pickle.load(open('G.pkl','rb'))
ASC=680; ERODE_R=11; DK=230/186.0; AK=None
FX={'b':.82,'d':.82,'p':.82,'q':.82,'g':.82,'m':.78,'k':.77,'v':.77,'x':.77,'z':.77,'f':.92,'j':.85,'l':.85,
    'n':.78,'u':.78,'c':.76,'t':.91,'i':.85}
SPEC_ASC=152/115.5*525
def ymap_inv(yp):
    # yp: target y_up before erosion -> source y_up
    half=ERODE_R; a=(525+2*ERODE_R)/525
    asc_k=(ASC+half-(525+half))/(SPEC_ASC-525)
    y=np.where(yp<-half, (yp+half)/DK,
      np.where(yp<=525+half, (yp+half)/a, 525+(yp-525-half)/asc_k))
    return y
def synth2(l):
    d=L[l]; fx=FX[l]
    x,y,w,h=d['box']; m=10
    S=XH_UNITS/115.5
    c=d['img'][y-m:y+h+m, x-m:x+w+m].astype(np.float32)
    big=cv2.resize(c,None,fx=S,fy=S,interpolation=cv2.INTER_CUBIC)
    big=cv2.GaussianBlur(big,(0,0),0.8*S)
    oy=(d['base']-(y-m))*S
    H,W=big.shape
    top=ASC+120; bot=-400
    newH=top-bot; newoy=top
    rows=np.arange(newH,dtype=np.float32)
    src_r=(oy-ymap_inv(newoy-rows)).astype(np.float32)
    newW=int(W*fx)+60
    src_c=((np.arange(newW,dtype=np.float32)-30)/fx).astype(np.float32)
    mx,my=np.meshgrid(src_c,src_r)
    out=cv2.remap(big,mx.astype(np.float32),my.astype(np.float32),cv2.INTER_LINEAR,borderValue=0)
    b=(out>128).astype(np.uint8)
    # isotropic thinning (done on an un-condensed copy is ideal; approximate with elliptical SE
    # whose horizontal radius is pre-shrunk by fx so that it is circular in source space)
    r=ERODE_R
    se=cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(2*max(1,int(round(r*fx)))+1,2*r+1))
    b=cv2.erode(b,se)
    dl=STEM_T-(STEM_S-2*r)*fx; k=int(round(abs(dl)))
    if k>=1:
        se=np.ones((1,k+1),np.uint8); b=cv2.dilate(b,se) if dl>0 else cv2.erode(b,se)
    if l in 'ij': b=fix_dot(b,newoy)
    return b,30,newoy
DOT=118
def fix_dot(b,oy):
    n,lab,st,cen=cv2.connectedComponentsWithStats(b)
    i=min(range(1,n),key=lambda i:st[i][1])
    cx=st[i][0]+st[i][2]/2; cy=st[i][1]+st[i][3]/2
    b[lab==i]=0
    # stem center for i
    if True:
        rows=b[int(oy-300)]; xs=np.nonzero(rows)[0]
        if len(xs): cx=(xs[0]+xs[-1])/2
    cy=oy-(ASC-DOT/2)  # dot top at ascender
    cv2.circle(b,(int(round(cx*4)),int(round(cy*4))),int(DOT/2*4),1,-1,lineType=cv2.LINE_AA,shift=2)
    return b
def traced(key,xh_px=None):
    b,ox,oy=G[key]
    return b,ox,oy
def sym_o(b,oy):
    ys,xs=np.nonzero(b); x0,x1=xs.min(),xs.max()
    cx=(x0+x1)/2; cy=oy-262.5
    H,W=b.shape
    yy,xx=np.mgrid[0:H,0:W]
    mx=np.where(xx<=cx,xx,np.round(2*cx-xx)).astype(int)
    my=np.where(yy>=cy,yy,np.round(2*cy-yy)).astype(int)
    ok=(mx>=0)&(mx<W)&(my>=0)&(my<H)
    out=np.zeros_like(b); out[ok]=b[my[ok],mx[ok]]
    return out

def bf_glyphs():
    r,res,t=sign_letters('bf')
    names=['B','a','y','s','h','o','r','e','F','w','y2']
    box=dict(zip(names,res[0]))
    out={}
    for nm in ['a','s','r','e','w']:
        g=box[nm]; b,ox,oy,S=raster(r,g,g[1]+g[3],g[3]); out[nm]=(b,ox,oy)
    g=box['h']; b,ox,oy,S=raster(r,g,g[1]+g[3],115); out['h']=(b,ox,oy)
    g=box['o']; b,ox,oy,S=raster(r,g,g[1]+g[3],115); out['o']=(sym_o(b,oy),ox,oy)
    g=box['y2']; b,ox,oy,S=raster(r,g,g[1]+114,114); out['y']=(b,ox,oy)
    return out
def smooth(b,sigma):
    f=cv2.GaussianBlur(b.astype(np.float32),(0,0),sigma)
    return (f>0.5).astype(np.uint8)
def all_glyphs():
    out={}
    for k,v in bf_glyphs().items():
        b,ox,oy=v; out[k]=(smooth(b,10),ox,oy,'sign')
    for l in 'bcdfgijklmnpqtuvxz':
        b,ox,oy=synth2(l); out[l]=(smooth(b,4),ox,oy,'synth')
    return out
