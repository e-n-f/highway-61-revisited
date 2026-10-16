import numpy as np,cv2,sys
from collections import deque
def silhouettes(path,open_k=4,thr=140):
    g=cv2.imread(path,0); H,W=g.shape
    ink=(g<thr).astype(np.uint8)
    body=ink.copy(); body[int(H*0.62):]=0                         # drop the table at the bottom of the page
    thick=cv2.morphologyEx(body,cv2.MORPH_OPEN,cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(open_k,open_k)))
    thick=cv2.morphologyEx(thick,cv2.MORPH_CLOSE,np.ones((5,5),np.uint8))   # bridge small gaps in the outline
    free=(1-thick).astype(np.uint8)
    n,lab=cv2.connectedComponents(free,connectivity=4)
    bg=lab[0,0]
    # adjacency between regions across the outline
    adj={i:set() for i in range(n)}
    k=cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(15,15))
    stats=cv2.connectedComponentsWithStats(free,connectivity=4)[2]
    for i in range(n):
        if stats[i][4]<300 or i==bg: continue
        x,y,w,h=stats[i][:4]; x0,y0=max(0,x-10),max(0,y-10); x1,y1=min(W,x+w+10),min(H,y+h+10)
        m=(lab[y0:y1,x0:x1]==i).astype(np.uint8); d=cv2.dilate(m,k)
        for j in np.unique(lab[y0:y1,x0:x1][(d>0)&(m==0)]):
            if j!=i: adj[i].add(j); adj[j].add(i)
    depth={bg:0}; q=deque([bg])
    while q:
        a=q.popleft()
        for b in adj[a]:
            if b not in depth and (stats[b][4]>=300): depth[b]=depth[a]+1; q.append(b)
    out=np.zeros((H,W),np.uint8)
    for r,dp in depth.items():
        if dp%2==1: out[lab==r]=1
    out=cv2.morphologyEx(out|thick&cv2.dilate(out,np.ones((9,9),np.uint8)),cv2.MORPH_CLOSE,np.ones((3,3),np.uint8))
    return out
if __name__=='__main__':
    p=sys.argv[1]; s=silhouettes(p)
    cv2.imwrite('sil_test.png',cv2.resize(255-s*255,None,fx=0.3,fy=0.3))
