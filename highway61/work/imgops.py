"""The few raster operations the scripts need, with numpy, scipy.ndimage and Pillow (they replace OpenCV calls, and
match them: the same structuring elements, border rules, connectivity and label order).

Images are numpy arrays, rows downwards; points are (x, y) in pixels."""
import numpy as np
from scipy import ndimage
from PIL import Image, ImageDraw, ImageFont


def imread_gray(fn):
    """8-bit greyscale image (cv2.imread(fn, 0))"""
    return np.array(Image.open(fn).convert('L'))


def imwrite(fn, im):
    """write an RGB array (channels in RGB order)"""
    Image.fromarray(np.asarray(im, np.uint8)).save(fn)


def imread_rgb(fn):
    return np.array(Image.open(fn).convert('RGB'))


def ellipse_kernel(w, h):
    """cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (w, h))"""
    k = np.zeros((h, w), bool); r, c = h // 2, w // 2
    for i in range(h):
        dy = i - r
        if abs(dy) <= r:
            dx = int(round(c * np.sqrt((r * r - dy * dy) / (r * r)))) if r else c
            k[i, max(c - dx, 0):min(c + dx + 1, w)] = True
    return k


def erode(a, k):
    """binary erosion; pixels outside the image count as set (cv2's default border for erosion)"""
    return ndimage.binary_erosion(np.asarray(a) > 0, structure=k, border_value=1)


def dilate(a, k):
    """binary dilation; pixels outside the image count as clear"""
    return ndimage.binary_dilation(np.asarray(a) > 0, structure=k, border_value=0)


def opening(a, k):
    """cv2.morphologyEx(a, cv2.MORPH_OPEN, k) for a 0/1 image (k with an odd size, anchored at its centre)"""
    return dilate(erode(a, k), k)


def components(a, connectivity=8):
    """cv2.connectedComponentsWithStats: (count including the background, labels, stats rows x y w h area)"""
    st = np.ones((3, 3), bool) if connectivity == 8 else ndimage.generate_binary_structure(2, 1)
    lab, n = ndimage.label(np.asarray(a) > 0, structure=st)
    if n:                       # numbered in the order OpenCV meets them: 2 x 2 blocks, row pairs from the top
        ys, xs = np.indices(lab.shape); key = (ys // 2) * lab.shape[1] + xs
        first = ndimage.minimum(key, lab, np.arange(1, n + 1))
        new = np.zeros(n + 1, int); new[1 + np.argsort(first, kind='stable')] = np.arange(1, n + 1)
        lab = new[lab]
    stats = np.zeros((n + 1, 5), int)
    stats[0] = (0, 0, lab.shape[1], lab.shape[0], int((lab == 0).sum()))
    area = np.bincount(lab.ravel(), minlength=n + 1)
    for i, sl in enumerate(ndimage.find_objects(lab), 1):
        stats[i] = (sl[1].start, sl[0].start, sl[1].stop - sl[1].start, sl[0].stop - sl[0].start, area[i])
    return n + 1, lab, stats


def boundary_points(a, min_points=0):
    """the pixels of a 0/1 image's contours as cv2.findContours(..., CHAIN_APPROX_NONE) gives them, as one unordered
    (x, y) array: set pixels with a clear 4-neighbour (or on the image's edge).  Contour pieces of min_points pixels
    or fewer (counted per 8-connected piece) are left out"""
    m = np.asarray(a) > 0
    b = m & ~ndimage.binary_erosion(m, structure=ndimage.generate_binary_structure(2, 1), border_value=0)
    if min_points:
        lab, n = ndimage.label(b, structure=np.ones((3, 3), bool))
        keep = np.bincount(lab.ravel(), minlength=n + 1) > min_points; keep[0] = False
        b = keep[lab]
    ys, xs = np.nonzero(b)
    return np.c_[xs, ys]


def gaussian(g, sigma):
    """cv2.GaussianBlur(g, (0, 0), sigma) for a float image: kernel size from sigma as OpenCV chooses it for
    floating-point images, reflected border (BORDER_REFLECT_101)"""
    ks = int(round(sigma * 8 + 1)) | 1
    return ndimage.gaussian_filter(np.asarray(g, float), sigma, mode='mirror', truncate=(ks // 2) / sigma)


def fill_poly(shape, polys):
    """0/1 mask of polygons (points in pixels, floats allowed), each filled with 1 in turn (as cv2.fillPoly called once
    per polygon; OpenCV also sets some pixels whose centres lie just outside an edge, under a pixel's width)"""
    im = Image.new('L', (shape[1], shape[0]), 0); d = ImageDraw.Draw(im)
    for p in polys:
        p = np.asarray(p, float)
        if len(p) > 2: d.polygon([tuple(q) for q in p], fill=1)
    return np.array(im)


SS = 4          # supersampling for anti-aliased drawing


def _mask(shape, draw):
    im = Image.new('L', (shape[1] * SS, shape[0] * SS), 0); draw(ImageDraw.Draw(im))
    return np.asarray(im.resize((shape[1], shape[0]), Image.BOX), float) / 255


def _blend(im, m, color):
    im[...] = np.round(im * (1 - m[..., None]) + np.array(color, float) * m[..., None]).astype(np.uint8)


def polylines(im, polys, color, closed=True, width=1):
    """anti-aliased lines on an RGB array (in place); points in pixels"""
    def draw(d):
        for p in polys:
            q = [tuple(x) for x in (np.asarray(p, float) + 0.5) * SS]
            if closed: q = q + q[:1]
            d.line(q, fill=255, width=width * SS)
    _blend(im, _mask(im.shape, draw), color)


def fill_polys(im, polys, color):
    """anti-aliased filled polygons (even-odd) on an RGB array, in place"""
    def draw(d):
        for p in polys: d.polygon([tuple(x) for x in (np.asarray(p, float) + 0.5) * SS], fill=255)
    _blend(im, _mask(im.shape, draw), color)


def line(im, p0, p1, color):
    """a one-pixel line between integer points (not anti-aliased), in place"""
    pi = Image.fromarray(im); ImageDraw.Draw(pi).line([tuple(p0), tuple(p1)], fill=tuple(color), width=1)
    im[...] = np.asarray(pi)


def text(im, xy, s, color=(0, 0, 0), size=14):
    """text with its baseline-left at xy, on an RGB array, in place"""
    try: f = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', size)
    except OSError: f = ImageFont.load_default()
    pi = Image.fromarray(im); ImageDraw.Draw(pi).text(xy, s, fill=tuple(color), font=f, anchor='ls')
    im[...] = np.asarray(pi)


def resize_area(a, s):
    """cv2.resize(a, None, fx=s, fy=s, interpolation=cv2.INTER_AREA) for an 8-bit single-channel image"""
    h, w = a.shape
    return np.asarray(Image.fromarray(np.asarray(a, np.uint8)).resize((int(round(w * s)), int(round(h * s))), Image.BOX))
