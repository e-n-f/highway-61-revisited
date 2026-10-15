"""specimen for Highway 61 Revisited <series>: the constructed capitals/figures over the 1966 font, and running text"""
import sys
from PIL import Image, ImageDraw, ImageFont, ImageChops
S = sys.argv[1] if len(sys.argv) > 1 else 'E'
F61 = '../outputs/Highway61Revisited-%s.otf' % S; F66 = '../outputs/HighwayEven1966-%s.ttf' % S
sz = 150
f61 = ImageFont.truetype(F61, sz); f66 = ImageFont.truetype(F66, sz); small = ImageFont.truetype(F61, 64)
lines = ['ABCDEFGHI', 'JKLMNOPQR', 'STUVWXYZ', '1234567890']
W, H = 1500, 4 * 190 + 420
img = Image.new('RGB', (W, H), 'white')
for i, t in enumerate(lines):
    y = 20 + i * 190
    a = Image.new('L', (W, 190), 0); ImageDraw.Draw(a).text((30, 0), t, font=f66, fill=255)
    b = Image.new('L', (W, 190), 0); ImageDraw.Draw(b).text((30, 0), t, font=f61, fill=255)
    # 1966 only: light blue, 1961 only: red, both: dark grey
    both = ImageChops.multiply(a, b)
    layer = Image.new('RGB', (W, 190), 'white')
    layer.paste((150, 190, 235), mask=a); layer.paste((220, 40, 40), mask=b); layer.paste((60, 60, 60), mask=both)
    img.paste(layer, (0, y))
d = ImageDraw.Draw(img)
d.text((30, 790), 'GOLDEN GATE BRIDGE  EXIT 432B', font=small, fill='black')
d.text((30, 880), 'San Francisco  Oakland Bay Bridge', font=small, fill='black')
d.text((30, 970), 'NORTH 101  Palo Alto 2  Junction 280', font=small, fill='black')
d.text((30, 1100), 'grey = both    red = 1961 construction only    blue = 1966 tracing only', font=ImageFont.truetype(F61, 36), fill=(90, 90, 90))
img.save('../outputs/specimen-1961-%s.png' % S)
