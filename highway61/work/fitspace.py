import json,string,subprocess,sys
from fontTools.ttLib import TTFont
s=sys.argv[1]
OT={'B':'SeriesB2000','C':'SeriesC2000','D':'SeriesD2000','E':'SeriesE2000','F':'SeriesF2000','EM':'SeriesModE2000'}[s]
OUT='/mnt/user-data/outputs/'+('Highway1966-EM.ttf' if s=='EM' else 'HighwayEven1966-%s.ttf'%s)
j=json.load(open('spacing_delta.json')); j[s]['lc']=0.0; json.dump(j,open('spacing_delta.json','w'))
subprocess.run(['python3','-W','ignore','buildfam.py',s],check=True,capture_output=True)
A=TTFont(OUT); B=TTFont('/mnt/user-data/outputs/fhwa2000/%s-Normal.otf'%OT)
cur=sum(A['hmtx'][c][0] for c in string.ascii_lowercase); tgt=sum(B['hmtx'][c][0] for c in string.ascii_lowercase)*690/285
d=round((tgt-cur)/52,1); j[s]['lc']=d; json.dump(j,open('spacing_delta.json','w'))
subprocess.run(['python3','-W','ignore','buildfam.py',s],check=True,capture_output=True)
print(s,'lc delta',d)
