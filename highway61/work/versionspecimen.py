"""The three versions of each series' lower case, one under another: python3 versionspecimen.py
-> ../../specimen/specimen-lc-versions.png

Main (lccondense.F_LC, spacing.LC_TIGHT: the a-z as wide as FHWA 2000's), Matched (one width factor and one spacing
factor per series giving FHWA's ratios to Series E Modified, s and w left out; lcmatch.py) and Natural (the
capitals' width factor, the spacing table scaled by the stroke); see variant.py."""
import os
from PIL import Image, ImageDraw, ImageFont
W = os.path.dirname(os.path.abspath(__file__)); DIST = os.path.join(W, '..', '..')
VERSIONS = [('', 'Main'), ('Matched', 'Matched'), ('Natural', 'Natural')]


def main():
    SER = ['A', 'B', 'C', 'D', 'E', 'F']; LINE = 150; H = 90 + len(VERSIONS) * 2 * LINE + 40
    im = Image.new('RGB', (3000, len(SER) * H + 20), 'white'); d = ImageDraw.Draw(im)
    lab = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', 36)
    small = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', 26)
    for i, S in enumerate(SER):
        y = i * H + 60
        d.text((20, y), 'Series ' + S, font=lab, fill=(80, 80, 80), anchor='ls')
        for j, (ps, name) in enumerate(VERSIONS):
            f = ImageFont.truetype(os.path.join(DIST, 'otf', 'Highway61Revisited%s-%s.otf' % (ps, S)), 140)
            yy = y + 30 + j * 2 * LINE
            d.text((20, yy + 70), name, font=small, fill=(150, 150, 150), anchor='ls')
            d.text((200, yy + 110), 'abcdefghijklmnopqrstuvwxyz', font=f, fill='black', anchor='ls')
            d.text((200, yy + 110 + LINE), 'Bayshore Fwy Business District', font=f, fill='black', anchor='ls')
    im.save(os.path.join(DIST, 'specimen', 'specimen-lc-versions.png'))


if __name__ == '__main__': main()
