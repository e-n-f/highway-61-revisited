import numpy as np,pickle
from scipy.optimize import least_squares
import geomc3 as g
G=pickle.load(open('geomfit.pkl','rb'))
def ell(p): return dict(a=p[2],b=p[3])
# candidate rings: (outer a,b), (inner a,b) and centre offset inner-outer (dx,dy)
C={}
for name,key in [('o','o'),('p bowl','p'),('b bowl','b'),('q bowl','q'),('d bowl','d')]:
    po,pi,oy,_=G[key]; C[name]=((po[2],po[3]),(pi[2],pi[3]),(pi[0]-po[0],-(pi[1]-po[1])))
def test(l,region,okind,ikind):
    Po=g.edge_pts(l,region[0],region[1],okind); Pi=g.edge_pts(l,region[0],region[1],ikind)
    out={}
    for name,((ao,bo),(ai,bi),(dx,dy)) in C.items():
        for mx in (1,-1):
            def res(p):
                cx,cy=p
                ro=np.sqrt(((Po[:,0]-cx)/ao)**2+((Po[:,1]-cy)/bo)**2)-1
                ri=np.sqrt(((Pi[:,0]-(cx+mx*dx))/ai)**2+((Pi[:,1]-(cy+dy))/bi)**2)-1
                return np.r_[ro*min(ao,bo),ri*min(ai,bi)]
            x0=[np.r_[Po[:,0],Pi[:,0]].mean(),np.r_[Po[:,1],Pi[:,1]].mean()]
            s=least_squares(res,x0)
            out[(name,mx)]=np.sqrt((s.fun**2).mean())
    best=sorted(out.items(),key=lambda kv:kv[1])[:3]
    return best,len(Po),len(Pi)
tests={
 'r arc':('r',((180,600),(380,600)),'top','bottom'),
 'n arch':('n',((180,500),(380,600)),'top','bottom'),
 'a top arc':('a',((60,470),(380,600)),'top','bottom'),
 'a bowl':('a',((0,330),(-30,330)),'bottom','top'),
 'f hook':('f',((120,330),(560,720)),'top','bottom'),
 't tail':('t',((140,400),(-30,150)),'bottom','top'),
 'j tail':('j',((0,180),(-300,-60)),'bottom','top'),
 'g tail':('g',((100,480),(-300,-40)),'bottom','top'),
 'u bowl':('u',((0,480),(-30,200)),'bottom','top'),
}
for k,(l,reg,ok,ik) in tests.items():
    best,no,ni=test(l,reg,ok,ik)
    print('%-10s'%k,no,ni,'  '.join('%s%s rms %.1f'%(n,'' if m==1 else '(mirror)',v) for (n,m),v in best))
print('--- refined regions')
tests2={
 'r arc':('r',((185,390),(380,600)),'top','bottom'),
 't tail':('t',((200,330),(-30,140)),'bottom','top'),
 'f hook':('f',((150,320),(600,720)),'top','bottom'),
 'g tail':('g',((120,400),(-300,-60)),'bottom','top'),
 'u bowl':('u',((120,400),(-30,150)),'bottom','top'),
 'a top arc':('a',((150,400),(420,600)),'top','bottom'),
}
for k,(l,reg,ok,ik) in tests2.items():
    best,no,ni=test(l,reg,ok,ik)
    print('%-10s'%k,no,ni,'  '.join('%s%s rms %.1f'%(n,'' if m==1 else '(mirror)',v) for (n,m),v in best))
