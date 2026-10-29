"""The standard spacing tables (1966 Standard Alphabets pp. 24-28, 34 and 36; 1983 pp. 29-33 and 41) against the fonts:
python3 spacingcheck.py [series...]

Capitals and figures: Table I (letter) / II (figure) give a code from the preceding character and the group of the
following one; Table VI gives the space for each code, the distance between the preceding character's extreme right
edge and the following one's extreme left edge.  The codes are the same in every series and both editions; the
spaces are the series stroke times 1.5, 1.2, 0.8, 0.4 (B C D E F, both editions), and for E(M) the Series E spaces
(1983: 1.3, 1.05, 0.7, 0.35 of its own stroke, height/5).  Series A has no table; it gets the same factors of its own stroke.

Lower case: the 1966 table (p. 36; E(M) only; the other series are given it scaled by spacing.lc_scale), in inches
for loop heights (= x-height) of 10, 12, 15 and 18 in, by preceding character class (capitals and lower case) and
following lower-case class.  The 1983 table is its 15-inch column in millimetres, but labelled as for a 375 mm loop
height, which made its spaces 1.6% too large for their height."""
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
# 1966 lower-case table (p. 36): inches, for loop heights of 10, 12, 15 and 18 in, by preceding character (capitals
# and lower case) and following lower-case class: acdegoq, bhiklmnpru, fw, j, st, vy, x, z (four sizes each).  The
# 1983 table is its 15-inch column in millimetres, labelled as for a 375 mm loop height (15 in is 381 mm).
LCROWS=['AWX','B','CEG','DOQR','F','HIMN','JU','KL','P','S','T','V','Y','Z','adghijlmnqu','bfkops','ce','r','tz','vy',
        'w','x']
LCCOLS=['acdegoq','bhiklmnpru','fw','j','st','vy','x','z']
LCSIZES=(10,12,15,18)
T1966={
'AWX':   '3.13 3.75 4.75 5.63 3.50 4.25 5.25 6.38 2.88 3.50 4.38 5.25 2.00 2.38 3.00 3.50 2.50 3.00 3.75 4.50 2.13 2.50 3.13 3.75 2.50 3.00 3.75 4.50 3.25 3.88 4.88 5.75',
'B':     '3.75 4.50 5.63 6.75 4.75 5.75 7.13 8.63 3.50 4.25 5.25 6.38 2.50 3.00 3.75 4.50 3.50 4.13 5.13 6.13 2.88 3.50 4.38 5.25 2.88 3.50 4.38 5.25 3.88 4.63 5.75 6.88',
'CEG':   '3.50 4.25 5.25 6.38 4.38 5.25 6.50 7.88 3.25 3.88 4.88 5.75 2.25 2.75 3.50 4.13 3.13 3.75 4.75 5.63 3.13 3.75 4.75 5.63 3.38 4.00 5.00 6.00 3.63 4.38 5.50 6.50',
'DOQR':  '3.50 4.13 5.13 6.13 4.50 5.38 6.75 8.00 3.38 4.00 5.00 6.00 2.25 2.63 3.25 3.88 3.13 3.75 4.75 5.63 3.13 3.75 4.75 5.63 3.38 4.00 5.00 6.00 3.50 4.25 5.25 6.38',
'F':     '1.63 2.00 2.50 3.00 3.38 4.00 5.00 6.00 2.13 2.50 3.13 3.75 1.63 2.00 2.50 3.00 1.88 2.25 2.75 3.38 1.88 2.25 2.75 3.38 2.00 2.38 3.00 3.50 2.50 3.00 3.75 4.50',
'HIMN':  '4.63 5.50 6.88 8.25 5.38 6.50 8.13 9.75 4.38 5.25 6.50 7.88 3.50 4.25 5.25 6.38 4.25 5.13 6.38 7.63 4.25 5.13 6.38 7.63 4.50 5.38 6.75 8.00 4.75 5.63 7.00 8.38',
'JU':    '4.50 5.38 6.75 8.00 4.75 5.75 7.13 8.63 4.25 5.13 6.38 7.63 3.50 4.13 5.13 6.13 4.00 4.75 6.00 7.13 4.00 4.75 6.00 7.13 4.13 5.00 6.25 7.50 4.63 5.50 6.88 8.25',
'KL':    '2.75 3.25 4.00 4.88 4.00 4.75 6.00 7.13 2.50 3.00 3.75 4.50 1.63 1.88 2.38 2.75 2.38 2.88 3.63 4.25 2.38 2.88 3.63 4.25 2.63 3.13 3.88 4.63 2.88 3.50 4.38 5.25',
'P':     '3.25 3.88 4.88 5.75 3.75 4.50 5.63 6.75 3.00 3.63 4.50 5.38 1.63 2.00 2.50 3.00 2.75 3.25 4.00 4.88 2.75 3.25 4.00 4.88 2.88 3.50 4.38 5.25 3.13 3.75 4.75 5.63',
'S':     '3.13 3.75 4.75 5.63 4.25 5.13 6.38 7.63 2.88 3.50 4.38 5.25 2.13 2.50 3.13 3.75 2.88 3.38 4.25 5.00 2.88 3.38 4.25 5.00 3.00 3.63 4.50 5.38 3.25 3.88 4.88 5.75',
'T':     '2.50 3.00 3.75 4.50 4.00 4.75 6.00 7.13 2.25 2.75 3.50 4.13 1.88 2.25 2.75 3.38 2.50 3.00 3.75 4.50 2.50 3.00 3.75 4.50 2.75 3.25 4.00 4.88 2.88 3.50 4.38 5.25',
'V':     '2.25 2.75 3.50 4.13 3.75 4.50 5.63 6.75 2.50 3.00 3.75 4.50 2.00 2.38 3.00 3.50 2.88 3.38 4.25 5.00 2.88 3.38 4.25 5.00 3.00 3.63 4.50 5.38 3.25 3.88 4.88 5.75',
'Y':     '1.63 2.00 2.50 3.00 3.38 4.00 5.00 6.00 2.13 2.50 3.13 3.75 1.50 1.75 2.25 2.63 1.63 1.88 2.38 2.88 1.63 1.88 2.38 2.88 1.88 2.25 2.75 3.38 2.75 3.25 4.00 4.88',
'Z':     '4.13 5.00 6.25 7.50 5.00 6.00 7.50 9.00 4.00 4.88 6.13 7.25 2.88 3.50 4.38 5.25 4.00 4.75 6.00 7.13 4.00 4.75 6.00 7.13 4.13 5.00 6.25 7.50 4.38 5.25 6.50 7.88',
'adghijlmnqu': '4.00 4.88 6.13 7.25 5.13 6.13 7.63 9.13 3.88 4.63 5.75 6.88 2.88 3.50 4.38 5.25 3.75 4.50 5.63 6.75 3.75 4.50 5.63 6.75 4.00 4.75 6.00 7.13 4.13 5.00 6.25 7.50',
'bfkops':'2.88 3.50 4.38 5.25 4.00 4.88 6.13 7.25 2.75 3.25 4.00 4.88 1.75 2.13 2.63 3.13 2.63 3.13 3.88 4.63 2.63 3.13 3.88 4.63 2.88 3.38 4.25 5.00 3.00 3.63 4.50 5.38',
'ce':    '3.13 3.75 4.75 5.63 4.13 5.00 6.25 7.50 2.88 3.50 4.38 5.25 2.00 2.38 3.00 3.50 2.88 3.38 4.25 5.00 2.88 3.38 4.25 5.00 3.00 3.63 4.50 5.38 3.25 3.88 4.88 5.75',
'r':     '2.00 2.38 3.00 3.50 3.13 3.75 4.75 5.63 1.88 2.25 2.75 3.38 0.67 1.00 1.25 1.50 1.63 2.00 2.50 3.00 1.63 2.00 2.50 3.00 1.88 2.25 2.75 3.38 2.13 2.50 3.13 3.75',
'tz':    '3.00 3.63 4.50 5.38 4.13 5.00 6.25 7.50 2.88 3.38 4.25 5.00 1.88 2.25 2.75 3.38 2.75 3.25 4.00 4.88 2.75 3.25 4.00 4.88 2.88 3.50 4.38 5.25 3.13 3.75 4.75 5.63',
'vy':    '2.63 3.13 3.88 4.63 3.75 4.50 5.63 6.75 2.38 2.88 3.63 4.25 1.50 1.75 2.25 2.63 2.25 2.75 3.50 4.13 2.25 2.75 3.50 4.13 2.50 3.00 3.75 4.50 2.75 3.25 4.00 4.88',
'w':     '2.75 3.25 4.00 4.88 3.75 4.50 5.63 6.75 2.50 3.00 3.75 4.50 1.63 1.88 2.38 2.75 2.38 2.88 3.63 4.25 2.38 2.88 3.63 4.25 2.63 3.13 3.88 4.63 2.88 3.38 4.25 5.00',
'x':     '2.88 3.38 4.25 5.00 3.88 4.63 5.75 6.88 2.63 3.13 3.88 4.63 1.63 2.00 2.50 3.00 2.50 3.00 3.75 4.50 2.50 3.00 3.75 4.50 2.75 3.25 4.00 4.88 2.88 3.50 4.38 5.25'}
def _lcratio():
    """space / loop height for each (row, class): the least-squares fit over the four sizes (the printed values are
    rounded to 1/8 in; rms residual 0.03 in).  The r -> j cell's 10-inch value (0.67) is left out: the others give
    0.83."""
    H=np.array(LCSIZES,float); out=[]
    for r in LCROWS:
        v=np.array([float(x) for x in T1966[r].split()]).reshape(8,4); row=[]
        for c in range(8):
            m=np.ones(4,bool)
            if r=='r' and LCCOLS[c]=='j': m[0]=False
            row.append(float(v[c][m]@H[m]/(H[m]@H[m])))
        out.append(row)
    return out
LCRATIO=_lcratio()
NAME={str(i):n for i,n in enumerate(['zero','one','two','three','four','five','six','seven','eight','nine'])}
def kerning(f):
    """{(left glyph, right glyph): x advance adjustment} from the GPOS pair lookups"""
    out={}
    if 'GPOS' not in f: return out
    for lk in f['GPOS'].table.LookupList.Lookup:
        for st in lk.SubTable:
            if lk.LookupType==9: st=st.ExtSubTable
            if getattr(st,'LookupType',lk.LookupType)!=2: continue
            if st.Format==1:
                for g,ps in zip(st.Coverage.glyphs,st.PairSet):
                    for r in ps.PairValueRecord: out.setdefault((g,r.SecondGlyph),getattr(r.Value1,'XAdvance',0) or 0)
            else:
                c1=st.ClassDef1.classDefs; c2=st.ClassDef2.classDefs
                for g in st.Coverage.glyphs:
                    rec=st.Class1Record[c1.get(g,0)]
                    for h,k in c2.items():
                        v=getattr(rec.Class2Record[k].Value1,'XAdvance',0) or 0
                        if v: out.setdefault((g,h),v)
    return out
def bearings(fn):
    """({char: (lsb, rsb)}, {(char, char): kern})"""
    f=TTFont(fn); gs=f.getGlyphSet(); hm=f['hmtx']; cm=f.getBestCmap(); out={}
    for ch in string.ascii_letters+string.digits:
        n=cm[ord(ch)]; p=BoundsPen(gs); gs[n].draw(p); b=p.bounds
        out[ch]=(b[0],hm[n][0]-b[2])
    kg=kerning(f); chars=string.ascii_letters+string.digits
    kern={(a,b):kg.get((cm[ord(a)],cm[ord(b)]),0) for a in chars for b in chars}
    return out,kern
def cap_pairs(S):
    """(preceding, following, standard gap, font gap) for every capital pair and figure pair, font units"""
    sb,kn=bearings(os.path.join(ROOT,'otf','Highway61Revisited-%s.otf'%S)); st=STROKE[S]; out=[]
    for a,codes in CODES.items():
        for b,g in GROUP.items(): out.append((a,b,FACTOR[codes[g]]*st,sb[a][1]+sb[b][0]+kn[(a,b)]))
    for a,codes in NCODES.items():
        for b,g in NGROUP.items(): out.append((a,b,FACTOR[codes[g]]*st,sb[a][1]+sb[b][0]+kn[(a,b)]))
    return out
def lc_pairs(S='EM'):
    """the same for the lower-case table; other series than E(M) scaled by spacing.lc_scale"""
    import spacing
    sc=spacing.lc_scale(S)
    sb,kn=bearings(os.path.join(ROOT,'otf','Highway61Revisited-%s.otf'%S)); out=[]
    for r,row in zip(LCROWS,LCRATIO):
        for cols,v in zip(LCCOLS,row):
            for a in r:
                for b in cols: out.append((a,b,v*XH*sc,sb[a][1]+sb[b][0]+kn[(a,b)]))
    return out
def report(pairs,label):
    d=np.array([f-s for _,_,s,f in pairs]); print('%-10s pairs %4d  font-standard: mean %+6.1f  mean|.| %5.1f  max|.| %5.1f'
          %(label,len(d),d.mean(),np.abs(d).mean(),np.abs(d).max()))
    return d
if __name__=='__main__':
    for S in (sys.argv[1:] or ['A','B','C','D','E','EM','F']):
        report(cap_pairs(S),S+' caps'); report(lc_pairs(S),S+' lc')
