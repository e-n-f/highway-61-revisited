"""The FHWA 2000 Standard Alphabets' embedded Type 1 fonts as OpenType:
python3 fhwa2otf.py [Alphabets.ps] -> ../fhwa2000/<name>.pfa (if the PostScript is given) and <name>.otf

The .pfa files are the %%BeginResource font sections of the publication's PostScript (pdftops), extracted verbatim
(see HANDOFF.md).  They are subsets: Series B..F 2000 and Series ModE 2000 have the whole alphabet, figures and
common punctuation; Highway D 98 and ME 98 only the letters the book sets in them.  Outlines and advance widths are
copied unchanged (CFF, 1000 units per em); the cmap comes from the glyph names."""
import glob, os
from fontTools.t1Lib import T1Font
from fontTools.fontBuilder import FontBuilder
from fontTools.pens.t2CharStringPen import T2CharStringPen
from fontTools.pens.boundsPen import BoundsPen
from fontTools.agl import toUnicode
D=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','fhwa2000')
def extract(ps):
    """write the Series*/Highway* font resources of the publication's PostScript as .pfa files"""
    import re
    s=open(ps,'rb').read()
    for m in re.finditer(rb'%%BeginResource: font (\S+)\r?\n(.*?)%%EndResource',s,re.S):
        base=m.group(1).decode().split('+')[-1]
        if base.startswith(('Series','Highway')):
            body=m.group(2); open(os.path.join(D,base+'.pfa'),'wb').write(body[body.find(b'%!'):])
def convert(fn):
    t=T1Font(fn); t.parse(); f=t.font; gs=t.getGlyphSet()
    name=f['FontName'].split('+')[-1]
    order=['.notdef']+sorted(n for n in f['CharStrings'].keys() if n!='.notdef')
    if 'space' not in order: order.append('space')
    cs,widths,cmap={},{},{}
    for n in order:
        pen=T2CharStringPen(None,None)
        if n in gs:
            gs[n].draw(pen); w=gs[n].width
        else: w=250                                    # space is missing from some subsets
        cs[n]=pen.getCharString(); widths[n]=w
        u=toUnicode(n)
        if len(u)==1: cmap[ord(u)]=n
    fb=FontBuilder(1000,isTTF=False); fb.setupGlyphOrder(order); fb.setupCharacterMap(cmap)
    fam=f['FontInfo'].get('FamilyName',name.split('-')[0])
    fb.setupCFF(name,{'FullName':name,'FamilyName':fam,'Weight':'Normal'},cs,{})
    bp={}
    for n in order:
        b=BoundsPen(None)
        if n in gs: gs[n].draw(b)
        bp[n]=b.bounds
    fb.setupHorizontalMetrics({n:(widths[n],(bp[n][0] if bp[n] else 0)) for n in order})
    fb.setupHorizontalHeader(ascent=920,descent=-80)
    fb.setupNameTable({'familyName':fam,'styleName':'Regular','psName':name,
                       'copyright':f['FontInfo'].get('Notice',''),'fullName':name})
    fb.setupOS2(sTypoAscender=920,sTypoDescender=-80,usWinAscent=920,usWinDescent=110)
    fb.setupPost(); fb.save(os.path.join(D,name+'.otf'))
    return name,len(order)
if __name__=='__main__':
    import sys
    if len(sys.argv)>1: extract(sys.argv[1])
    for fn in sorted(glob.glob(os.path.join(D,'*.pfa'))):
        print(*convert(fn))
