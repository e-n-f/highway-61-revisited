import numpy as np,yaml,json,pickle,cv2
import geomc2,geomc3,geomc4,geomc5,geomc6,geomc7,geomc8,geomc9,narrowmk,assemble
from geomc3 import master,runs_at,sc,term_line,TM
R=lambda v: float(round(float(v),1))
def x0of(l): b,oy=master(l); return int(np.nonzero(b.any(0))[0][0])
def ell(cx,cy,a,b,l): return {'cx':R(cx-x0of(l)),'cy':R(cy),'rx':R(a),'ry':R(b)}
def ell_rows(p,oy,l): return ell(p[0],oy-p[1],p[2],p[3],l)          # fit_super params (x,row,a,b)
def stem(l,run): return {'left':R(run[0]-x0of(l)),'right':R(run[1]-x0of(l)),'centre':R(sc(run)-x0of(l))}
def term(l,which):
    t=term_line(l,which)
    if t is None: return None
    ang=np.degrees(np.arctan2(t['u'][1],t['u'][0]))
    if ang<=-90: ang+=180
    if ang>90: ang-=180
    return {'x':R(t['p'][0]-x0of(l)),'y':R(t['p'][1]),'angle_deg':R(ang)}
def profile(key,n=24):
    pr=geomc6.PROF[key] if key in geomc6.PROF else geomc5.PROF[key]
    th=pr['th']; out=[]
    for k in range(n):
        a=2*np.pi*k/n; i=np.argmin(np.abs(th-a)); out.append([R(np.degrees(a)),R(pr['t'][i])])
    return out
L={}
# o: brush ring
po,pi,oy,_=geomc2.G['o']
L['o']={'construction':'even brush ring (constant-width pen around an elliptical centre line)',
        'outer_ellipse':ell_rows(po,oy,'o'),'inner_ellipse_1966':ell_rows(pi,oy,'o')}
# bowls b d p q
for l in 'bdpq':
    po,pi,oy=geomc5.BOWL[l]; st=runs_at(l,300)[0] if l in 'bp' else runs_at(l,300)[-1]
    up=l in 'bd'
    L[l]={'construction':'bowl = inner ellipse + outward pen (thickness profile) ; stem rectangle with slanted %s cut'%('ascender' if up else 'descender'),
          'bowl_outer_ellipse':ell_rows(po,oy,l),'bowl_inner_ellipse':ell_rows(pi,oy,l),
          'pen_profile_deg_thickness':profile(l),'stem':stem(l,st),
          'stem_cut':term(l,'top' if up else 'bottom'),'stem_square_end_y':0 if up else 525}
# a
po,pi,oy=geomc5.BOWL['a']; s=runs_at('a',150)[-1]
c=lambda e: ell(e[0],e[1],e[2],e[3],'a')
L['a']={'construction':'bowl = inner ellipse + pen profile; top arc = inside TOP, outside BIG and outside TOP_INNER; stem on the right',
        'bowl_outer_ellipse':ell_rows(po,oy,'a'),'bowl_inner_ellipse':ell_rows(pi,oy,'a'),'bowl_pen_profile':profile('a'),
        'arc_top_ellipse':c(geomc9.A_TOP),'arc_big_ellipse':c(geomc9.A_BIG),'arc_top_inner_ellipse':c(geomc9.A_TIN),
        'stem':stem('a',s),'arc_terminal':term('a','top'),
        'note':'arc ellipses come from your screenshot markup (converted to font units)'}
# c
po,pi,oy,_=geomc6.FIT['c']
L['c']={'construction':'inner ellipse + pen profile; opening between the two terminal cuts',
        'outer_ellipse':ell_rows(po,oy,'c'),'inner_ellipse':ell_rows(pi,oy,'c'),'pen_profile':profile('c'),
        'terminals':[ {'x':R(Lc['c'][0]-x0of('c')),'y':R(oy-Lc['c'][1]),'angle_deg':R(np.degrees(np.arctan2(-Lc['u'][1],Lc['u'][0])))}
                      for Lc in geomc3.ML['c'] if 3<np.degrees(np.arctan2(abs(Lc['u'][1]),abs(Lc['u'][0])))<87]}
# e
L['e']={'construction':'top half and bottom half from two different ellipse pairs (inner + pen profile); crossbar; aperture cut',
        'top_outer_ellipse':ell_rows(geomc7.po_t,geomc7.oy,'e'),'top_inner_ellipse':ell_rows(geomc7.pi_t,geomc7.oy,'e'),
        'bottom_outer_ellipse':ell_rows(geomc7.po_b,geomc7.oy,'e'),'bottom_inner_ellipse':ell_rows(geomc7.pi_b,geomc7.oy,'e'),
        'top_pen_profile':profile('e_top'),'bottom_pen_profile':profile('e_bot'),
        'crossbar_bottom_y':R(geomc7.YB),'crossbar_top_y':R(geomc7.YT),
        'aperture_terminal':term('e','bottom')}
# g
po,pi,oy=geomc5.BOWL['q']
L['g']={'construction':'bowl = the q bowl; stem; tail = inner ellipse (0.81x outer, touching it at upper left) + pen profile; terminal cut',
        'bowl':'same as q, moved to the g stem','stem':stem('g',runs_at('g',300)[-1]),
        'tail_outer_ellipse':ell_rows(geomc8.po,geomc8.oy,'g'),'tail_inner_ellipse':ell_rows(geomc8.pi,geomc8.oy,'g'),
        'tail_inner_scale':[R(geomc8.ei[0]*100)/100,R(geomc8.ei[1]*100)/100],'tail_pen_profile':profile('g_tail'),
        'tail_terminal':term('g','bottom')}
# n h m r u
po,pi,oy,_=geomc6.FIT['n']
arch={'outer_ellipse':ell_rows(po,oy,'n'),'inner_ellipse':ell_rows(pi,oy,'n'),'pen_profile':profile('n')}
L['n']={'construction':'left stem; arch = inner ellipse spanning the counter + pen profile (kept above the inner ellipse centre); right leg',
        'stems':[stem('n',r_) for r_ in runs_at('n',150)],'arch':arch}
po,pi,oy,_=geomc6.FIT['h']
L['h']={'construction':'as n, with ascender stem','stems':[stem('h',r_) for r_ in runs_at('h',150)],
        'arch':{'outer_ellipse':ell_rows(po,oy,'h'),'inner_ellipse':ell_rows(pi,oy,'h'),'pen_profile':profile('h')},
        'ascender_cut':term('h','top')}
L['m']={'construction':'two copies of the n arch, one per pair of stems','stems':[stem('m',r_) for r_ in runs_at('m',150)],'arch':'same as n'}
L['r']={'construction':'stem + the n arch (counter width = n stem spacing), ended by the terminal cut',
        'stem':stem('r',runs_at('r',150)[0]),'arch':'same as n','terminal':term('r','right')}
po,pi,oy,_=geomc6.FIT['u']
L['u']={'construction':'n arch construction upside down','stems':[stem('u',r_) for r_ in runs_at('u',450)],
        'bowl':{'outer_ellipse':ell_rows(po,oy,'u'),'inner_ellipse':ell_rows(pi,oy,'u'),'pen_profile':profile('u')}}
# t f
c=lambda e,l: ell(e[0],e[1],e[2],e[3],l)
L['t']={'construction':'stem (slanted top cut) + crossbar; tail inside BLUE, outside ORANGE; RED rounds the stem lower-left',
        'stem':stem('t',max(runs_at('t',250),key=lambda r_:r_[1]-r_[0])),'crossbar':{'x0':R(36-x0of('t')),'x1':R(385-x0of('t')),'top_y':527,'thickness_at_full_weight':131},
        'red_ellipse':c(geomc9.T_RED,'t'),'blue_ellipse':c(geomc9.T_BLUE,'t'),'orange_ellipse':c(geomc9.T_ORG,'t'),
        'ascender_cut':term('t','top'),'tail_terminal':term('t','bottom')}
L['f']={'construction':'stem + crossbar; hook inside OUTER, outside INNER, ending at the crossbar end',
        'stem':stem('f',max(runs_at('f',250),key=lambda r_:r_[1]-r_[0])),'crossbar':{'x0':R(36-x0of('f')),'x1':R(314-x0of('f')),'top_y':525,'thickness_at_full_weight':141},
        'outer_ellipse':c(geomc9.F_OUT,'f'),'inner_ellipse':c(geomc9.F_IN,'f')}
# j: p-bowl ring part
rp=geomc4.rp_of('p'); pst=sc(runs_at('p',300)[0])
L['j']={'construction':'stem + lower-right part of the p bowl ring (two ellipses), vertical cut at the left; round dot',
        'stem':stem('j',max(runs_at('j',300),key=lambda r_:r_[1]-r_[0])),
        'ring':'p bowl (outer/inner ellipses of p), right extreme on the j stem centre, bottom on the j bottom','offset_fitted':[R(v) for v in geomc4.OFF.get('j',(0,0))]}
# straight letters: record stroke centre lines by drawing at full width
rec={}
class Spy(geomc3.Canvas):
    def diag(s,p0,p1,y0,y1):
        rec.setdefault(cur,[]).append({'p0':[R(p0[0]-x0of(cur)),R(p0[1])],'p1':[R(p1[0]-x0of(cur)),R(p1[1])],'y_from':R(y0),'y_to':R(y1)})
        return super().diag(p0,p1,y0,y1)
orig=geomc3.Canvas; geomc3.Canvas=Spy
for cur in 'kvwxy': geomc3.build(cur,1.0,140)
geomc3.Canvas=orig
for l in 'kvwxy':
    seen=[];segs=[]
    for d in rec[l]:
        if (d['p0'],d['p1']) not in seen: seen.append((d['p0'],d['p1'])); segs.append(d)
    L[l]={'construction':'straight strokes of width T along these centre lines (direction p0->p1), cut flat at the y limits',
          'strokes':segs}
L['k']['stem']=stem('k',runs_at('k',150)[0]); L['k']['ascender_cut']=term('k','top')
L['w']['join_trim']={'amount_fraction_of_T':0.16,'over_height':230,'note':'outer edge pulled in at the bottoms and centre top'}
L['i']={'construction':'stem 0..525 + round dot','stem':stem('i',runs_at('i',150)[0])}
L['l']={'construction':'stem with slanted ascender cut','stem':stem('l',runs_at('l',150)[0]),'ascender_cut':term('l','top')}
L['z']={'construction':'two bars + diagonal (perpendicular width T) starting 0.25 bar-thickness into each bar; bars cut flush with the diagonal',
        'notes':'bar thickness = 1966 bar thickness x T/140; x extent from the 1966 z'}
L['s']={'construction':'traced from the 1966 drawing and narrowed (not constructed)'}
series={}
import json as _j
sd=_j.load(open('spacing_delta.json'))
for s,(f,T) in list(narrowmk.SER.items())+[('E',(1.0,123)),('E(M)',(1.0,140))]:
    key={'E(M)':'EMC'}.get(s,s)
    series[s]={'width_factor_f':f,'stroke_T':T,'spacing_delta_per_side':sd.get(key)}
doc={'about':{'units':'font units (1000/em); x measured from the left edge of the 1966 letter, y measured up from the baseline',
              'metrics':{'x_height':525,'cap_height':690,'ascender':700,'descender_approx':-260},
              'full_weight_stroke_TM':140,
              'how_widths_are_made':'x positions of centre lines are multiplied by f; strokes are drawn at width T; heights stay as drawn; pen profiles are scaled by T/140 (or so the stroke joining a stem equals T)',
              'pen_profile_format':'[angle around the inner ellipse in degrees (0 = right, 90 = top), thickness out to the outer edge]'},
     'series':series,'letters':{k:L[k] for k in sorted(L)}}
yaml.safe_dump(doc,open('/mnt/user-data/outputs/lowercase-skeleton.yaml','w'),sort_keys=False,width=140,allow_unicode=True)
print('written', sum(len(str(v)) for v in L.values()))
