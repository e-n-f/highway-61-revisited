import pickle,numpy as np
from build import smooth
SER={'B':85,'C':95,'D':102,'E':113,'F':123}
OUTER={'b':'L','h':'L','k':'L','l':'L','t':'L','f':'L','d':'R'}
def runs(row):
    xs=np.nonzero(row)[0]
    if len(xs)==0: return []
    sp=np.where(np.diff(xs)!=1)[0]
    return list(zip(np.r_[xs[0],xs[sp+1]],np.r_[xs[sp],xs[-1]]))
for s,T in SER.items():
    d=pickle.load(open('lcE_%s.pkl'%s,'rb'))
    for l,side in OUTER.items():
        b,ox,oy,tag=d[l]; b=b.copy()
        # reference ascender run at y=600
        ra=runs(b[int(oy-600)])
        if not ra: continue
        a0,a1=min(ra,key=lambda r:abs((r[0]+r[1])/2-np.mean([(x+y)/2 for x,y in ra])))
        cs=[]
        for y in range(180,400,4):
            for x0,x1 in runs(b[int(oy-y)]):
                if x1>=a0-40 and x0<=a1+40 and (x1-x0+1)<=T+15: cs.append((x0+x1)/2)
        if not cs: print('skip',s,l); continue
        c=float(np.median(cs)); L0=int(round(c-T/2)); R0=int(round(c+T/2))
        top=max(y for y in range(500,760) if b[int(oy-y)].any())
        if l=='f':
            for y in range(540,top):
                rr=[r for r in runs(b[int(oy-y)]) if min(r[1],R0)-max(r[0],L0)>=0.5*T]
                if not rr or (rr[0][1]-rr[0][0]+1)>T+15: top=y-12; break
        for y in range(4 if l=='l' else 300,top+1):
            r=int(round(oy-y)); row=b[r]
            for x0,x1 in runs(row):
                if min(x1,R0)-max(x0,L0) < 0.5*T: continue
                w=x1-x0+1
                if w<=T+40:
                    row[x0:x1+1]=0
                    if w>=T-6 or y<400: row[L0:R0]=1
                    else:  # slanted top cut: keep width, align to the side it hugged
                        if (x0-L0)>(R0-x1): row[R0-w:R0]=1
                        else: row[L0:L0+w]=1
                else:   # merged with arch/bowl/crossbar: fix outer edge only
                    if l=='t' or l=='f': continue
                    if side=='L' and x0<L0+40: row[x0:x1+1]=0; row[L0:x1+1]=1
                    if side=='R' and x1>R0-40: row[x0:x1+1]=0; row[x0:R0]=1
        d[l]=(smooth(b,2),ox,oy,tag)
    pickle.dump(d,open('lcE_%s.pkl'%s,'wb'))
    b,ox,oy,_=d['h']
    print(s,[ (lambda r:(r[0][0],r[0][1]))(runs(b[int(oy-y)])) for y in (200,400,560,640)])
