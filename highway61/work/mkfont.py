import pickle,numpy as np,cv2
from vec import *
from fontTools.fontBuilder import FontBuilder
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.cu2qu import curve_to_quadratic
A=pickle.load(open('A.pkl','rb'))
A.update(pickle.load(open('A2.pkl','rb')))
DIG=dict(zip('0123456789',['zero','one','two','three','four','five','six','seven','eight','nine'])); DIG.update({'-':'hyphen','.':'period',',':'comma','!':'exclam',"'":'quotesingle','"':'quotedbl','\u2018':'quoteleft','\u2019':'quoteright','\u201c':'quotedblleft','\u201d':'quotedblright'})
S,R,D=85,58,36
SB={ # lsb, rsb
 'a':(R,S),'b':(S,R),'c':(R,R-8),'d':(R,S),'e':(R,R),'f':(48,30),'g':(R,S),'h':(S,S),'i':(S,S),
 'j':(30,S),'k':(S,D),'l':(S,S),'m':(S,S),'n':(S,S),'o':(R,R),'p':(S,R),'q':(R,S),'r':(S,40),
 's':(R,R),'t':(45,48),'u':(S,S),'v':(D,D),'w':(D,D),'x':(D,D),'y':(D,D),'z':(48,48),
 'A':(20,20),'B':(S,R),'C':(R,50),'D':(S,R),'E':(S,45),'F':(S,40),'G':(R,70),'H':(S,S),'I':(S,S),
 'J':(40,S),'K':(S,22),'L':(S,30),'M':(S,S),'N':(S,S),'O':(R,R),'P':(S,30),'Q':(R,R),'R':(S,30),
 'S':(R,38),'T':(25,25),'U':(S,S),'V':(20,20),'W':(20,20),'X':(20,20),'Y':(20,20),'Z':(45,45),
 '-':(40,40),'.':(45,45),',':(45,45),'!':(70,70),"'":(40,40),'"':(40,40),'\u2018':(40,40),'\u2019':(40,40),'\u201c':(40,40),'\u201d':(40,40),'0':(R,R),'1':(60,S),'2':(50,50),'3':(50,R),'4':(35,55),'5':(R,R),'6':(R,R),'7':(45,35),'8':(R,R),'9':(R,R)}
YS=[0,525,680,690,-230]
FORCE_FLAT=False
EXTRA_SB={}
SMW=7
FORCE_SET={'A':[690]}
def smooth_closed(P,w=2):
    out=P.copy()
    for k in range(1,w+1): out=out+np.roll(P,k,0)+np.roll(P,-k,0)
    return out/(2*w+1)
def process(C):
    P=resample(C,3.0)
    n=len(P); cs=corners(smooth_closed(P,3))
    if FORCE_FLAT:
        for yy in FORCE_FLAT:
            on=np.abs(P[:,1]-yy)<2.0
            if not on.any() or on.all(): continue
            idx=np.nonzero(on)[0]
            # run boundaries (cyclic)
            for i in idx:
                if not on[(i-1)%n] or not on[(i+1)%n]:
                    cs=[c for c in cs if min(abs(i-c),n-abs(i-c))>6]+[int(i)]
        cs=sorted(set(cs))
    # snap corner coordinates
    for j in range(len(cs)):
        a,b=cs[j],cs[(j+1)%len(cs)]
        pa,pb=P[a],P[b]; d=pb-pa; L=np.hypot(*d)
        if L<30: continue
        seg=P[np.arange(a,(b if b>a else b+n)+1)%n]
        dev=np.abs(np.cross(d,seg-pa))/L
        if dev.max()>2.5: continue
        if abs(d[0])<0.05*abs(d[1]):
            m=(pa[0]+pb[0])/2; P[a,0]=P[b,0]=m
        elif abs(d[1])<0.05*abs(d[0]):
            m=(pa[1]+pb[1])/2; P[a,1]=P[b,1]=m
    # corner-preserving heavy smoothing inside each segment
    if len(cs)>=1:
        Q=P.copy(); W=SMW
        for j in range(len(cs)):
            a=cs[j]; b=cs[(j+1)%len(cs)]
            if b<=a: b+=n
            idx=np.arange(a,b+1)%n
            seg=P[idx]
            for t in range(1,len(seg)-1):
                w=min(W,t,len(seg)-1-t)
                Q[idx[t]]=seg[t-w:t+w+1].mean(0)
        P=Q
    else:
        Q=P.copy()
        for k in range(1,7): Q=Q+np.roll(P,k,0)+np.roll(P,-k,0)
        P=Q/13
    for c in cs:
        y=min(YS,key=lambda v:abs(v-P[c,1]))
        if abs(P[c,1]-y)<9: P[c,1]=y
    return postsnap(contour_segments(P))
def setpt(segs,i,which,new):
    # move endpoint of seg i ('s' start or 'e' end), keep neighbours consistent
    n=len(segs)
    def mv(j,end,new):
        s=list(segs[j])
        if end=='s':
            d=np.asarray(new)-s[1]; s[1]=np.asarray(new,float)
            if s[0]=='C': s[2]=s[2]+d
        else:
            if s[0]=='L': d=np.asarray(new)-s[2]; s[2]=np.asarray(new,float)
            else: d=np.asarray(new)-s[4]; s[4]=np.asarray(new,float); s[3]=s[3]+d
        segs[j]=tuple(s)
    if which=='s': mv(i,'s',new); mv((i-1)%n,'e',new)
    else: mv(i,'e',new); mv((i+1)%n,'s',new)
def endp(s): return s[2] if s[0]=='L' else s[4]
def postsnap(segs):
    segs=list(segs); n=len(segs)
    for it in range(2):
        for i in range(n):
            s=segs[i]
            if s[0]!='L': continue
            a,b=s[1],s[2]; d=b-a; L=np.hypot(*d)
            if L<20: continue
            if abs(d[0])<0.06*abs(d[1]):
                # include collinear vertical neighbours
                xs=[(a[0]+b[0])/2]; w=[L]; grp=[i]
                j=(i+1)%n
                while segs[j][0]=='L' and j!=i:
                    dd=segs[j][2]-segs[j][1]
                    if abs(dd[0])<0.06*abs(dd[1])+1 and abs(segs[j][1][0]-b[0])<4: xs.append((segs[j][1][0]+segs[j][2][0])/2); w.append(np.hypot(*dd)); grp.append(j); j=(j+1)%n
                    else: break
                x=np.average(xs,weights=w)
                for g in grp:
                    sg=segs[g]; setpt(segs,g,'s',(x,sg[1][1])); sg=segs[g]; setpt(segs,g,'e',(x,sg[2][1]))
            elif abs(d[1])<0.06*abs(d[0]):
                y=(a[1]+b[1])/2
                yy=min(YS,key=lambda v:abs(v-y))
                if abs(y-yy)<10: y=yy
                setpt(segs,i,'s',(a[0],y)); setpt(segs,i,'e',(b[0],y))
    # merge consecutive collinear lines
    out=[]
    for s in segs:
        if out and s[0]=='L' and out[-1][0]=='L':
            a=out[-1][1]; b=s[2]; m=s[1]; d=b-a; L=np.hypot(*d)
            if L>0 and abs(d[0]*(m-a)[1]-d[1]*(m-a)[0])/L<0.6:
                out[-1]=('L',a,b); continue
        out.append(s)
    return out
def glyph_outlines(b,oy):
    cnts,hier=cv2.findContours(b.copy(),cv2.RETR_CCOMP,cv2.CHAIN_APPROX_NONE)
    res=[]
    for i,c in enumerate(cnts):
        if len(c)<20 or cv2.contourArea(c)<250: continue
        pts=c[:,0,:].astype(float)
        P=np.c_[pts[:,0],oy-pts[:,1]]
        area=0.5*np.sum(P[:,0]*np.roll(P[:,1],-1)-np.roll(P[:,0],-1)*P[:,1])
        hole=hier[0][i][3]!=-1
        # outer: clockwise (area<0); hole: ccw
        if (not hole and area>0) or (hole and area<0): P=P[::-1]
        res.append(process(P))
    return res
def draw(pen,segs,dx):
    def T(p): return (int(round(p[0]+dx)),int(round(p[1])))
    pen.moveTo(T(segs[0][1]))
    for s in segs:
        if s[0]=='L': pen.lineTo(T(s[2]))
        else:
            q=curve_to_quadratic([tuple(p) for p in s[1:]],1.0)
            pts=[T(p) for p in q[1:]]
            # drop degenerate
            pen.qCurveTo(*pts)
    pen.closePath()
def build(path,family="Highway Gothic 1966 Narrow",style="Regular",src=None,sbscale=1.0,desc=None,lcadd=0,lcscale=None):
    lcs=sbscale if lcscale is None else lcscale
    G=src or A
    glyphs={}; metrics={}
    CH='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.,!\'"\u2018\u2019\u201c\u201d-'
    gn=lambda c: DIG.get(c,c)
    order=['.notdef','space']+[gn(c) for c in CH]
    for g in ['.notdef','space']:
        pen=TTGlyphPen(None)
        if g=='.notdef':
            for rect in [((50,0),(50,680),(400,680),(400,0)),((100,50),(350,50),(350,630),(100,630))]:
                pen.moveTo(rect[0]); [pen.lineTo(p) for p in rect[1:]]; pen.closePath()
            metrics[g]=(450,50)
        else: metrics[g]=(int(300*sbscale),0)
        glyphs[g]=pen.glyph()
    for l in CH:
        b,ox,oy,_s=G[l]
        global FORCE_FLAT
        FORCE_FLAT=FORCE_SET.get(l,False)
        outs=glyph_outlines(b,oy)
        FORCE_FLAT=False
        xs=[p[0] for segs in outs for s in segs for p in s[1:]]
        _c=np.nonzero(b.any(0))[0]; xmin=float(_c[0]); xmax=float(_c[-1]+1)
        lsb,rsb=[v*(lcs if l.islower() else sbscale)+(lcadd if l.islower() else 0)+EXTRA_SB.get(l,0) for v in SB[l]]; dx=lsb-xmin
        pen=TTGlyphPen(None)
        for segs in outs: draw(pen,segs,dx)
        glyphs[gn(l)]=pen.glyph()
        adv=int(round(xmax-xmin+lsb+rsb))
        metrics[gn(l)]=(adv,lsb)
    fb=FontBuilder(1000,isTTF=True)
    fb.setupGlyphOrder(order)
    cmap={ord(c):gn(c) for c in CH}; cmap[32]='space'
    fb.setupCharacterMap(cmap)
    fb.setupGlyf(glyphs)
    # recompute lsb from glyf bbox
    gs=fb.font['glyf']
    hm={}
    for g in order:
        gl=gs[g]; gl.recalcBounds(gs)
        hm[g]=(metrics[g][0], getattr(gl,'xMin',0) if gl.numberOfContours else 0)
    fb.setupHorizontalMetrics(hm)
    fb.setupHorizontalHeader(ascent=900,descent=-340)
    fb.setupNameTable({"familyName":family,"styleName":style,
        "uniqueFontIdentifier":family.replace(' ','')+"-"+style,
        "fullName":family+" "+style,"psName":family.replace(' ','')+"-"+style,
        "version":"Version 1.000",
        "copyright":"Outlines traced from FHWA Standard Alphabets for Highway Signs (1966) lower case sheets and photographs of two surviving narrow-lowercase signs in Palo Alto, CA.",
        "description":desc or "1960s narrow lowercase highway alphabet with 1966 Series D capitals and figures. Lowercase a e h o r s w y traced from signs; other lowercase proportionally narrowed and lightened from the 1966 standard-width lower case."})
    fb.setupOS2(sTypoAscender=900,sTypoDescender=-340,sTypoLineGap=0,usWinAscent=900,usWinDescent=340,
                sxHeight=525,sCapHeight=690,achVendID="NONE",fsType=0)
    fb.setupPost(); fb.setupDummyDSIG() if False else None
    fb.save(path)
if __name__=='__main__':
    build('/home/claude/w/test.ttf')
    print('built')
