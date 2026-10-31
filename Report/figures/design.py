"""Regenerate fig_design.png from the original task screens in the AIRIPA slides."""
from pathlib import Path
from io import BytesIO
import pymupdf
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

folder = Path(__file__).resolve().parent
slides = pymupdf.open(folder.parents[1] / "Slide/7_Pagan_AIRIPA2026.pdf")
fig = plt.figure(figsize=(10, 7.2), facecolor="white")
ax = fig.add_axes([0, 0, 1, 1])
ax.set(xlim=(0, 10), ylim=(0, 7.2))
ax.axis("off")

# Session sequences follow the preregistration; screenshots are unchanged.
rows = [
    (5.8, "A  T0", ["16 arithmetic\nproblems", "State worry\nand perceived\ncompetence", "16 arithmetic\nproblems", "Trait math anxiety\nand academic\nbuoyancy"]),
    (4.2, "B  T1", ["16 problems\nwith choice", "Break", "16 problems\nwith choice", "Trait math avoidance\nand self-competence"]),
]
for y, label, boxes in rows:
    ax.text(.25, y + 1.03, label, fontsize=13, weight="bold")
    for j, text in enumerate(boxes):
        x = .25 + j * 2.45
        ax.add_patch(Rectangle((x, y), 2.15, .84, facecolor="#f1f3f4", edgecolor="#555555", linewidth=.8))
        ax.text(x + 1.075, y + .42, text, ha="center", va="center", fontsize=10)
        if j < 3:
            ax.annotate("", (x + 2.44, y + .42), (x + 2.17, y + .42), arrowprops=dict(arrowstyle="->", lw=.9))

for pos, page, title in [(0, 4, "C  Choice screen"), (1, 5, "D  Hint screen")]:
    # Embedded-image order differs from the visual order on the slides.
    xref = slides[page].get_images()[2 if page == 4 else 1][0]
    screen = Image.open(BytesIO(slides.extract_image(xref)["image"]))
    panel = fig.add_axes([.025 + pos * .5, .03, .45, .47])
    panel.imshow(screen)
    panel.axis("off")
    panel.set_title(title, loc="left", fontsize=13, weight="bold", pad=10)

fig.savefig(folder / "fig_design.png", dpi=250, facecolor="white")
plt.close(fig)
