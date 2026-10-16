"""The standard spacing tables (1966 Standard Alphabets pp. 24-28 and 34; 1983 pp. 29-33 and 41) against the fonts:
python3 spacingcheck.py [series...]

Capitals and figures: Table I (letter) / II (figure) give a code from the preceding character and the group of the
following one; Table VI gives the space for each code, the distance between the preceding character's extreme right
edge and the following one's extreme left edge.  The codes are the same in every series and both editions; the
spaces are the series stroke times 1.5, 1.2, 0.8, 0.4 (B C D E F, both editions), and for E(M) the Series E spaces
(1983: 1.3, 1.05, 0.7, 0.35 of its own stroke, height/5).  Series A has no table; it is given Series B's rule.

Lower case (E(M) only): the 1983 table, in millimetres at a 375 mm loop height (= x-height), by preceding character
class (capitals and lower case) and following lower-case class.  The 1966 table agrees with it (inches)."""
import os, sys, string, numpy as np
from fontTools.ttLib import TTFont
from fontTools.pens.boundsPen import BoundsPen
W=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.join(W,'..','..')
CAP=690.0; XH=517.5
# Table I: preceding letter -> codes for following groups (B D E F H I K L M N P R U | C G O Q S X Z | A J T V W Y)
CODES=dict(A=(2,2,4),B=(1,2,2),C=(2,2,3),D=(1,2,2),E=(2,2,3),F=(2,2,3),G=(1,2,2),H=(1,1,2),I=(1,1,2),J=(1,1,2),
           K=(2,2,3),L=(2,2,4),M=(1,1,2),N=(1,1,2),O=(1,2,2),P=(1,2,2),Q=(1,2,2),R=(1,2,2),S=(1,2,2),T=(2,2,4),
           U=(1,1,2),V=(2,2,4),W=(2,2,4),X=(2,2,3),Y=(2,2,4),Z=(2,2,3))
GROUP={c:0 for c in 'BDEFHIKLMNPRU'}; GROUP.update({c:1 for c in 'CGOQSXZ'}); GROUP.update({c:2 for c in 'AJTVWY'})
# Table II: preceding figure -> codes for following groups (1 5 | 2 3 6 8 9 0 | 4 7)
NCODES={'1':(1,1,2),'2':(1,2,2),'3':(1,2,2),'4':(2,2,4),'5':(1,2,2),'6':(1,2,2),'7':(2,2,4),'8':(1,2,2),'9':(1,2,2),
        '0':(1,2,2)}
NGROUP={c:0 for c in '15'}; NGROUP.update({c:1 for c in '236890'}); NGROUP.update({c:2 for c in '47'})
FACTOR={1:1.5,2:1.2,3:0.8,4:0.4}
STROKE={'A':64.69,'B':86.25,'C':97.03,'D':107.81,'E':118.59,'F':129.38,'EM':118.59}   # EM: Series E's spaces
# 1983 lower-case table at 375 mm loop height; columns: acdegoq, bhiklmnpru, fw, j, st, vy, x, z
LCROWS=['AWX','B','CEG','DOQR','F','HIMN','JU','KL','P','S','T','V','Y','Z','adghijlmnqu','bfkops','ce','r','tz','vy',
        'w','x']
LCCOLS=['acdegoq','bhiklmnpru','fw','j','st','vy','x','z']
LC375=[[121,133,111,76,95,80,95,124],[143,181,133,95,130,111,111,146],[133,165,124,89,121,121,127,140],
       [130,171,127,82,121,121,127,133],[64,127,80,64,70,70,76,95],[175,206,165,133,162,162,171,178],
       [171,181,162,130,152,152,159,175],[102,152,95,60,92,92,98,111],[124,143,114,64,102,102,111,121],
       [121,162,111,80,108,108,114,124],[95,152,89,70,95,95,102,111],[89,143,95,76,108,108,114,124],
       [64,127,80,57,60,60,70,102],[159,190,156,111,152,152,159,165],[156,194,146,111,143,143,152,159],
       [111,156,102,67,98,98,108,114],[121,159,111,76,108,108,114,124],[76,121,70,32,64,64,70,80],
       [114,159,108,70,102,102,111,121],[98,143,92,57,89,89,95,102],[102,143,95,60,92,92,98,108],
       [108,146,98,64,95,95,102,111]]
NAME={str(i):n for i,n in enumerate(['zero','one','two','three','four','five','six','seven','eight','nine'])}
def bearings(fn):
    f=TTFont(fn); gs=f.getGlyphSet(); hm=f['hmtx']; cm=f.getBestCmap(); out={}
    for ch in string.ascii_letters+string.digits:
        n=cm[ord(ch)]; p=BoundsPen(gs); gs[n].draw(p); b=p.bounds
        out[ch]=(b[0],hm[n][0]-b[2])
    return out
def cap_pairs(S):
    """(preceding, following, standard gap, font gap) for every capital pair and figure pair, font units"""
    sb=bearings(os.path.join(ROOT,'Highway61Revisited-%s.otf'%S)); st=STROKE[S]; out=[]
    for a,codes in CODES.items():
        for b,g in GROUP.items(): out.append((a,b,FACTOR[codes[g]]*st,sb[a][1]+sb[b][0]))
    for a,codes in NCODES.items():
        for b,g in NGROUP.items(): out.append((a,b,FACTOR[codes[g]]*st,sb[a][1]+sb[b][0]))
    return out
def lc_pairs(S='EM'):
    sb=bearings(os.path.join(ROOT,'Highway61Revisited-%s.otf'%S)); out=[]
    for r,row in zip(LCROWS,LC375):
        for cols,v in zip(LCCOLS,row):
            for a in r:
                for b in cols: out.append((a,b,v*XH/375,sb[a][1]+sb[b][0]))
    return out
def report(pairs,label):
    d=np.array([f-s for _,_,s,f in pairs]); print('%-10s pairs %4d  font-standard: mean %+6.1f  mean|.| %5.1f  max|.| %5.1f'
          %(label,len(d),d.mean(),np.abs(d).mean(),np.abs(d).max()))
    return d
if __name__=='__main__':
    for S in (sys.argv[1:] or ['A','B','C','D','E','EM','F']):
        report(cap_pairs(S),S+' caps')
    report(lc_pairs('EM'),'EM lc')
