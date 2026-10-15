import numpy as np, dim1961 as p5, dimE1 as e1, dimE2 as e2, dimE3 as e3, dimE4 as e4, dimE6 as e6
def poly(fn):
    r=fn(); polys=r[0] if isinstance(r,tuple) else r
    return p5.raster(polys)
REG={'A':e1.letter_A,'B':e1.letter_B,'C':lambda:e1.letter_C()[0],'D':e1.letter_D,'E':e1.letter_E,'F':e1.letter_F,
     'G':lambda:e2.letter_G()[0],'H':e2.letter_H,'I':e2.letter_I,'J':lambda:e2.letter_J()[0],'K':lambda:e2.letter_K()[0],'L':e2.letter_L,
     'M':e3.letter_M,'N':e3.letter_N,'O':e3.letter_O,'P':e3.letter_P,'Q':e3.letter_Q,'R':lambda:e3.letter_R()[0],
     'S':lambda:e4.letter_S()[0],'T':e4.letter_T,'U':lambda:e4.letter_U()[0],'V':e4.letter_V,'W':e4.letter_W,'X':lambda:e4.letter_X()[0],
     'Y':lambda:poly(p5.Y),'Z':lambda:p5.raster(p5.Z()[0]),'1':lambda:p5.raster(p5.one()[0]),'2':lambda:p5.two()[0],'3':lambda:p5.three()[0],
     '4':lambda:poly(p5.four),'5':lambda:e6.fig5()[0],'6':lambda:e6.fig6()[0],'7':lambda:e6.fig7()[0],'8':lambda:e6.fig8()[0],
     '9':lambda:e6.fig9()[0],'0':e6.fig0}
