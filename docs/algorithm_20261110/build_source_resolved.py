"""Occurrence-aware schematic; SVG source, official VectorCraft PDF/native export.

Short synthetic bases illustrate an inversion plus a flanking substitution.
REF core AGTCCCGA has reverse complement TCGGGACT; TT changes to AC.
These drawing-scale sequences are below production evidence thresholds and are
not an executable benchmark or a literal representation of the oriP locus.
"""
from pathlib import Path
from html import escape
from math import atan2, sin, cos

OUT = Path(__file__).resolve().parent
K, G, LIGHT = '#272727', '#737373', '#EAEAEA'
B, R = '#008280', '#631879'
s = ['<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="58.5mm" viewBox="0 0 1200 390">',
     '<title>Occurrence-aware allele extraction from graph paths</title>']

def rect(x,y,w,h,c,stroke='none',sw=1):
    s.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{c}" stroke="{stroke}" stroke-width="{sw}"/>')

def path(d,c=K,w=1.5,fill='none',dash=''):
    s.append(f'<path d="{d}" fill="{fill}" stroke="{c}" stroke-width="{w}" stroke-linecap="round" stroke-linejoin="round"'+(f' stroke-dasharray="{dash}"' if dash else '')+'/>')

def line(x,y,xx,yy,c=K,w=1.5,dash=''):
    path(f'M{x},{y} L{xx},{yy}',c,w,dash=dash)

def arrow(x,y,xx,yy,c=K,w=1.5):
    a=atan2(yy-y,xx-x);bx=xx-7*cos(a);by=yy-7*sin(a)
    line(x,y,xx,yy,c,w)
    path(f'M{xx},{yy} L{bx-3*sin(a)},{by+3*cos(a)} M{xx},{yy} L{bx+3*sin(a)},{by-3*cos(a)}',c,w)

def text(x,y,t,size=14,c=K,bold=False,anchor='start'):
    s.append(f'<text x="{x}" y="{y}" font-family="Helvetica" font-size="{size}" font-weight="{700 if bold else 400}" fill="{c}" text-anchor="{anchor}">{escape(t)}</text>')

def node(x,y,t,c=B,w=32,h=23):
    rect(x,y,w,h,c);text(x+w/2,y+16.5,t,15,'white',anchor='middle')

def seqbox(x,y,w,t,col=K):
    rect(x,y,w,25,'white',col,1.2)
    text(x+w/2,y+18,t,14,K,anchor='middle')

rect(0,0,1200,390,'white')
text(16,30,'a',25,K,True)
text(52,30,'Occurrence-aware allele extraction',22,K,True)
line(895,24,916,24,B,2);text(925,29,'Reference',14,G)
line(1031,24,1052,24,R,2);text(1061,29,'Sample',14,G)
text(52,82,'Graph paths',18,K,True)

text(385,82,'Mapping + occurrences',18,K,True)
text(803,82,'Complete alleles',18,K,True)

# Rectangular base labels and curved links inherit only the supplied visual grammar.
# The synthetic graph stores the two alleles as forward branches; sequence RC is recognized downstream.
seqbox(48,180,48,'GGA')
seqbox(112,137,29,'TT',B);seqbox(157,137,119,'AGTCCCGA',B)
seqbox(112,223,29,'AC',R);seqbox(157,223,119,'TCGGGACT',R)
seqbox(300,180,48,'CTA')
path('M96,192 C105,192 101,149 112,149',B,1.4)
line(141,149,157,149,B,1.4)
path('M276,149 C290,149 286,192 300,192',B,1.4)
path('M96,192 C105,192 101,235 112,235',R,1.4)
line(141,235,157,235,R,1.4)
path('M276,235 C290,235 286,192 300,192',R,1.4)
arrow(174,125,258,125,B,1.2)
arrow(174,260,258,260,R,1.2)

arrow(353,149,373,149,K,1.3)
path('M353,235 L365,235 L365,294',K,1.2)
arrow(365,294,382,294,K,1.2)

# Mapping uses full sample paths; this drawing abbreviates the surrounding context.
text(385,122,'minimap2',14,K)
text(385,154,'Reference',13,G);text(385,180,'Sample',13,G)
for yy,col in [(149,B),(175,R)]:
    for a,b in [(459,473),(498,535),(560,612),(637,656)]:line(a,yy,b,yy,col,1.2)
    for xx in [473,535,612]:rect(xx,yy-5,25,10,col)
for xx in [485,547,624]:line(xx,155,xx,169,G,0.8)
text(385,237,'Node visits',14,K)
# A shared synthetic node GGA at different reference/query path positions.
seqbox(414,269,48,'GGA',B);text(438,319,'r',13,G,anchor='middle')
seqbox(571,309,48,'GGA',R);text(595,359,'q',13,G,anchor='middle')
line(466,282,567,322,G,1,dash='3 4')
arrow(418,258,458,258,B,1.1);arrow(575,298,615,298,R,1.1)
path('M675,162 L711,162 L711,230',K,1.3)
path('M675,321 L711,321 L711,230',K,1.3)
s.append(f'<circle cx="711" cy="230" r="3.5" fill="{K}"/>')
arrow(711,230,755,230,K,1.3)

# Black flanks abstract eligible surrounding anchors; the short displayed bases
# illustrate sequence relationships, not actual 31-bp production anchors.
text(770,160,'REF',13,G);text(770,240,'ALT',13,G)
for yy in [143,223]:
    rect(813,yy,10,25,K);rect(1153,yy,10,25,K)
    line(823,yy+12,834,yy+12,G,1)
    line(1136,yy+12,1153,yy+12,G,1)
path('M813,168 L813,223 L823,223 L823,168 Z','none',0,LIGHT)
path('M1153,168 L1153,223 L1163,223 L1163,168 Z','none',0,LIGHT)
seqbox(834,143,26,'TT',B);seqbox(834,223,26,'AC',R)
seqbox(873,143,122,'AGTC',B);seqbox(1010,143,126,'CCGA',B)
seqbox(873,223,122,'TCGG',R);seqbox(1010,223,126,'GACT',R)
for yy,col in [(155,B),(235,R)]:
    line(860,yy,873,yy,col,1.1);line(995,yy,1010,yy,col,1.1)
# The two four-base halves reverse their order and orientation.
line(934,170,1073,221,G,0.9,dash='3 4')
line(1073,170,934,221,G,0.9,dash='3 4')
arrow(883,132,985,132,B,1.1);arrow(1020,132,1126,132,B,1.1)
arrow(985,260,883,260,R,1.1);arrow(1126,260,1020,260,R,1.1)
rect(824,281,9,9,K);text(845,291,'Anchors',13,G)
arrow(851,300,851,318,K,1.1)
text(851,343,'Allele calls',14,K,True,'middle')
arrow(1083,273,1083,296,K,1.1)
text(1083,319,'Classify + project',14,K,True,'middle')
text(1083,343,'INV_COMPLEX',15,R,True,'middle')
text(1083,373,'bcftools norm · replay',13,G,anchor='middle')
s.append('</svg>')
(OUT/'PathAllele_panel_a.svg').write_text('\n'.join(s)+'\n')
print(OUT/'PathAllele_panel_a.svg')
