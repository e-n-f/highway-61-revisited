from src import *
FX={'b':.82,'d':.82,'p':.82,'q':.82,'g':.82,'m':.78,'k':.77,'v':.77,'x':.77,'z':.77,'f':.92,'j':.85,'l':.85,
    # for comparison only
    'n':.78,'h':.78,'u':.78,'a':.82,'o':.87,'e':.87,'c':.76,'w':.76,'y':.78,'r':.95,'t':.91,'s':.73,'i':.85}
def sign_glyphs():
    G={}
    r,res,t=sign_letters('bf')
    names=['B','a','y','s','h','o','r','e','F','w','y2']
    for nm,g in zip(names,res[0]):
        xh=115
        base=g[1]+g[3] if nm not in ('y','y2') else g[1]+xh
        if nm in ('B','F'): base=g[1]+g[3]
        b,ox,oy,S=raster(r,g,base,xh)
        G['bf_'+nm]=(b,ox,oy)
    r,res,t=sign_letters('bd')
    for line,word in zip(res,['Business','District']):
        xh=float(np.median([g[3] for g,c in zip(line,word) if c in 'usnerc']))
        for k,(nm,g) in enumerate(zip(word,line)):
            b,ox,oy,S=raster(r,g,g[1]+g[3],xh)
            G['bd%d_%s%d'%(word=='District',nm,k)]=(b,ox,oy)
    return G
def sym_o(b):
    ys,xs=np.nonzero(b); x0,x1,y0,y1=xs.min(),xs.max(),ys.min(),ys.max()
    cx=(x0+x1)/2; cy=(y0+y1)/2
    out=b.copy()
    H,W=b.shape
    for yy in range(H):
        for xx in range(W):
            # map to bottom-left quadrant
            mx=xx if xx<=cx else int(round(2*cx-xx))
            my=yy if yy>=cy else int(round(2*cy-yy))
            out[yy,xx]=b[my,mx] if 0<=mx<W and 0<=my<H else 0
    return out
