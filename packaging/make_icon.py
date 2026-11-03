"""Draw the SPOT icon (1024 px PNG): two homologous chromosomes with a pair of orthologous loci joined."""

import sys
from pathlib import Path

from PIL import Image, ImageDraw

S = 1024
out = Path(sys.argv[1])
img = Image.new("RGBA", (S, S), (0, 0, 0, 0))
d = ImageDraw.Draw(img)
# rounded square background (macOS icon grid: 824 px body, 100 px margin)
m, r = 100, 185
d.rounded_rectangle([m, m, S - m, S - m], r, fill=(23, 42, 69))
# chromosomes
bar_h, x0, x1 = 64, 210, 814
for y, color in ((360, (229, 234, 240)), (664, (160, 178, 199))):
    d.rounded_rectangle([x0, y - bar_h // 2, x1, y + bar_h // 2], bar_h // 2, fill=color)
# syntenic link (band) between the two loci
link = Image.new("RGBA", (S, S), (0, 0, 0, 0))
ld = ImageDraw.Draw(link)
ld.polygon([(470, 392), (566, 392), (640, 632), (544, 632)], fill=(193, 18, 31, 120))
img = Image.alpha_composite(img, link)
d = ImageDraw.Draw(img)
# loci
for cx, cy in ((518, 360), (592, 664)):
    d.ellipse([cx - 58, cy - 58, cx + 58, cy + 58], fill=(193, 18, 31), outline=(255, 255, 255), width=12)
out.parent.mkdir(parents=True, exist_ok=True)
img.save(out)
print(out)
