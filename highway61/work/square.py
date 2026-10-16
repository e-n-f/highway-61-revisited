import pickle,numpy as np
def runs(row):
    xs=np.nonzero(row)[0]
    if len(xs)==0: return []
    sp=np.where(np.diff(xs)!=1)[0]
    return list(zip(np.r_[xs[0],xs[sp+1]],np.r_[xs[sp],xs[-1]]))
def R(oy,y): return int(round(oy-y))
def col_extent(b,oy,x):
    ys=np.nonzero(b[:,x])[0]; return oy-ys.max(), oy-ys.min()   # (lowest y, highest y)
def square_stem(b,oy,yref,pick,end):
    rr=runs(b[R(oy,yref)]); L,Rt=rr[0] if pick=='first' else rr[-1]
    c=(L+Rt)//2
    # find the stem's end along its centre column, near the requested end
    col=b[:,c]; ys=oy-np.nonzero(col)[0]
    if end=='bottom':
        seg=ys[ys<yref]; yend=seg.min(); rows=range(R(oy,yend+25),R(oy,yend)+1)
    else:
        seg=ys[(ys>yref)&(ys<yref+400)]; yend=seg.max(); rows=range(R(oy,yend),R(oy,yend-25)+1)
    for r in rows: b[r,L:Rt+1]=1
    return L,Rt,yend
def square_bar(b,oy,stem_y=200,probe=None):
    rr=runs(b[R(oy,stem_y)]); L,Rt=rr[0]
    # bar vertical extent from a column just right of the stem
    x=Rt+25; lo,hi=col_extent(b,oy,x)
    col=b[:,x]; r=R(oy,500)
    assert col[r]
    a=r
    while col[a-1]: a-=1
    z=r
    while col[z+1]: z+=1
    yt,yb=oy-a,oy-z
    mid=R(oy,(yb+yt)/2); r0=runs(b[mid]); xl=min(a for a,_ in r0); xr=max(z for _,z in r0)
    for r in range(R(oy,yt),R(oy,yb)+1): b[r,xl:xr+1]=1
    # clear fillets in the inner corners (outside the stem, 15 units above/below the bar)
    for y in list(range(int(yt)+1,int(yt)+16))+list(range(int(yb)-15,int(yb))):
        r=R(oy,y)
        b[r,xl:L]=0; b[r,Rt+1:xr+1]=0
    return (yb,yt,xl,xr)
def extend_edge(b,oy,arm,ytop,fit=(380,470)):
    ys=np.arange(*fit); xs=[]
    for y in ys:
        rr=runs(b[R(oy,y)]); a,z=rr[0] if arm=='left' else rr[-1]
        xs.append(z if arm=='left' else a)
    k=np.polyfit(ys,xs,1)
    for y in range(fit[1],int(ytop)+1):
        r=R(oy,y); rr=runs(b[r])
        if not rr: continue
        x=int(round(np.polyval(k,y)))
        if arm=='left':
            a,z=rr[0]
            if x>z: b[r,z:x+1]=1
        else:
            a,z=rr[-1]
            if x<a: b[r,x:a+1]=1
def apply(d):
    for l,pick,yr in [('a','last',230),('b','first',260),('d','last',260),('u','last',260)]:
        b,ox,oy,t=d[l]; b=b.copy(); print(l,square_stem(b,oy,yr,pick,'bottom')); d[l]=(b,ox,oy,t)
    for l,pick in [('g','last'),('m','first'),('n','first'),('p','first'),('q','last'),('r','first')]:
        b,ox,oy,t=d[l]; b=b.copy(); print(l,square_stem(b,oy,300,pick,'top')); d[l]=(b,ox,oy,t)
    for l in 'ft':
        b,ox,oy,t=d[l]; b=b.copy(); print(l,square_bar(b,oy)); d[l]=(b,ox,oy,t)
    b,ox,oy,t=d['x']; b=b.copy()
    ytop=oy-np.nonzero(b.any(1))[0].min()
    extend_edge(b,oy,'left',ytop); extend_edge(b,oy,'right',ytop); d['x']=(b,ox,oy,t)
    return d
if __name__=='__main__':
    d=pickle.load(open('hr_lc_clean.pkl','rb'))
    d=apply(d)
    pickle.dump(d,open('hr_lc_sq.pkl','wb'))
