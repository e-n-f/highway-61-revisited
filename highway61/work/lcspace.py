"""Choose a series' lower-case spacing factor (spacing.LC_TIGHT) so that its alphabet (sum of the a-z advances) is
as wide as the FHWA 2000 font's, both in cap heights:  python3 lcspace.py <series...>

The advances are the build's (build1961.py imported for its glyph boxes and starting bearings, nothing written)
with the side bearings spacing.solve gives for a trial lower-case scale; the total is linear in the scale, so two
trials and a secant step (repeated once for the rounding) settle it.  Needs the FHWA fonts (fhwa2otf.py)."""
import sys, string, importlib
import spacing, lcwidth, lccondense
GEO=None
def total(S,cur,box,lc):
    sb,_=spacing.solve(S,cur,lc,geo=GEO)
    return sum(round(sb[c][0]+box[c][1]-box[c][0]+sb[c][1]) for c in string.ascii_lowercase)/690
if __name__=='__main__':
    for S in sys.argv[1:]:
        sys.argv=['build1961.py',S]; sys.modules.pop('build1961',None); b=importlib.import_module('build1961')
        GEO=b._geo
        T=lcwidth.target(S); base=lccondense.stroke(S)/lccondense.V_EM
        t0=total(S,b._cur,b._box,base); t1=total(S,b._cur,b._box,base*0.5)
        k=1+(T-t0)*(0.5-1)/(t1-t0); t2=total(S,b._cur,b._box,base*k)
        k=k+(T-t2)*(0.5-1)/(t1-t0); k=round(k*1000)/1000                  # once more for the rounding
        print(S,'target %.3f  table-scaled %.3f  LC_TIGHT %.3f -> %.3f'%(T,t0,k,total(S,b._cur,b._box,base*k)),flush=True)
