import numpy as np,pickle
M=pickle.load(open('hr_lc_sq.pkl','rb'))
def top_at(b,oy,x):
    ys=np.nonzero(b[:,x])[0]; return oy-ys.min()
def stem(b,oy,y=600):
    xs=np.nonzero(b[int(round(oy-y))])[0]; return xs[0],xs[-1]
def cut_params(l):
    b,ox,oy,_=M[l]; L,R=stem(b,oy); return top_at(b,oy,L+3),top_at(b,oy,R-3),R-L
def fix(b,oy,l):
    yl_m,yr_m,w_m=cut_params(l); dy=yr_m-yl_m
    L,R=stem(b,oy)
    ytop=max(top_at(b,oy,x) for x in range(L,R+1))
    hi,lo=(ytop,ytop-dy) if dy>0 else (ytop+dy,ytop)  # right/left
    for x in range(L,R+1):
        t=(x-L)/max(1,R-L)
        yc=(ytop-dy)+dy*t if dy>0 else ytop+dy*t   # cut height at column x
        for y in range(int(min(hi,lo))-30, int(ytop)+15):
            r=int(round(oy-y))
            b[r,x]=1 if y<=yc else 0
    return b
if __name__=='__main__':
    for l in 'bdhkl': print(l,cut_params(l))
    d=pickle.load(open('lcE_B.pkl','rb'))
    for l in 'bdhkl':
        b,ox,oy,t=d[l]; d[l]=(fix(b.copy(),oy,l),ox,oy,t)
    pickle.dump(d,open('lcE_B.pkl','wb'))
