"""report near-miss tangencies: joins between outline segments whose directions differ by 0.05-12 degrees"""
import numpy as np, vgeom, outline
def tan(s, end):
    if s[0] == 'L': d = s[2] - s[1]
    else:
        p = s[2] if end else s[1]; r = p - s[3]; d = np.array([-r[1], r[0]]) * np.sign(s[5])
    return d / np.linalg.norm(d)
if __name__ == '__main__':
    V = vgeom.install()
    for ch, g in V.items():
        for segs in outline.contours(g):
            for i in range(len(segs)):
                a, b = segs[i - 1], segs[i]
                t = np.degrees(np.arccos(np.clip(tan(a, 1) @ tan(b, 0), -1, 1)))
                if a[0] == 'A' and b[0] == 'A' and t < 12:
                    # the circles' own tangency gap: |r1 -/+ r2| - |c1 c2| (internal if the arcs turn the same way)
                    d = np.linalg.norm(a[3] - b[3])
                    gap = (abs(a[4] - b[4]) if np.sign(a[5]) == np.sign(b[5]) else a[4] + b[4]) - d
                    if abs(gap) > 0.05:
                        print(ch, 'arcs r=%.1f/%.1f at (%.1f, %.1f): gap %.2f units, %.2f deg' % (a[4], b[4], *b[1], gap, t))
                elif 0.05 < t < 12:
                    print(ch, '%s->%s at (%.1f, %.1f): %.2f deg' % (a[0], b[0], *b[1], t))
