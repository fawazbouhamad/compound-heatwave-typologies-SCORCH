"""Save complete figures and their unlettered panels.

Panel boxes use fractions of the full image, measured from the upper left.
Panel groups follow the manuscript layout. Shared headings and legends are
repeated where needed, without resampling the panel pixels. Figure C1's
colour scale is drawn at each standalone panel's height.
"""
from __future__ import annotations

import io
import os
from pathlib import Path
import re
import shutil
import subprocess
from uuid import uuid4

from PIL import Image, ImageChops

Image.MAX_IMAGE_PIXELS = None


def _panel(box, **parts):
    return {"box": box, **parts}


def _grid(xs, ys):
    return {
        chr(97 + row * (len(xs) - 1) + col):
        _panel((xs[col], ys[row], xs[col + 1], ys[row + 1]))
        for row in range(len(ys) - 1) for col in range(len(xs) - 1)
    }


PANEL_LAYOUTS = {
    "figure_1": {},
    "figure_2": _grid([0, 0.5, 1], [0, 1]),
    "figure_3": {
        "a": _panel((0, 0, 0.5, 0.337)),
        "b": _panel((0.5, 0, 1, 0.337)),
        "c": _panel((0, 0.337, 0.5, 0.671)),
        "d": _panel((0.5, 0.337, 1, 0.671)),
        "e": _panel((0.24, 0.671, 0.78, 1)),
    },
    # Unlettered types and snapshots use left-to-right, top-to-bottom order.
    "figure_4": {
        "a": _panel((0, 0, 0.247, 1)),
        "b": _panel((0.253, 0, 0.497, 1)),
        "c": _panel((0.503, 0, 0.747, 1)),
        "d": _panel((0.753, 0, 1, 1)),
    },
    "figure_5": {
        "a": _panel((0, 0, 0.515, 0.895), footer=(0, 0.895, 1, 1)),
        "b": _panel((0.515, 0, 1, 0.895), footer=(0, 0.895, 1, 1)),
    },
    # Each letter in Figure 6 identifies an event containing two daily maps.
    "figure_6": {
        "a": _panel((0, 0, 1, 0.522), footer=(0, 0.934, 1, 1)),
        "b": _panel((0, 0.522, 1, 0.934), header=(0, 0, 1, 0.13),
                    footer=(0, 0.934, 1, 1), keep_width=True),
    },
    "figure_7": {
        **{letter: _panel((x0, y0, x1, y1), header=(0, 0, 1, 0.10),
                          footer=(0, 0.93, 1, 1))
           for letter, x0, x1, y0, y1 in [
               ("a", 0, 0.345, 0.10, 0.545),
               ("b", 0.345, 0.673, 0.10, 0.545),
               ("c", 0.673, 1, 0.10, 0.545),
               ("d", 0, 0.345, 0.545, 0.93),
               ("e", 0.345, 0.673, 0.545, 0.93),
               ("f", 0.673, 1, 0.545, 0.93),
           ]},
    },
    # The manuscript's a/b/c panels are entire rows of four heatwave types.
    "figure_8": {
        "a": _panel((0, 0, 1, 0.335)),
        "b": _panel((0, 0.335, 1, 0.675), header=(0, 0, 1, 0.028),
                    keep_width=True),
        "c": _panel((0, 0.675, 1, 1), header=(0, 0, 1, 0.028),
                    keep_width=True),
    },
    "figure_9": {
        "a": _panel((0, 0, 0.333, 0.90), footer=(0, 0.90, 1, 1)),
        # The final x tick in (b) extends beneath (c)'s y-axis label.
        # Preserve each label and exclude only its neighbour's margin text.
        "b": _panel((0.333, 0, 0.67, 0.90), footer=(0, 0.90, 1, 1),
                    exclude=[(0.6585, 0, 0.67, 0.77)]),
        "c": _panel((0.6585, 0, 1, 0.90), footer=(0, 0.90, 1, 1),
                    exclude=[(0.6585, 0.77, 0.67, 0.90)]),
    },
    "figure_10": _grid([0, 1], [0, 0.46, 1]),
    "figure_11": _grid([0, 0.5, 1], [0, 0.5, 1]),
    "figure_12": {},
    "figure_A1": _grid([0, 0.53, 1], [0, 0.333, 0.644, 1]),
    "figure_A2": {
        letter: _panel((x0, y0, x1, y1), footer=(0, 0.935, 1, 1))
        for letter, x0, x1, y0, y1 in [
            ("a", 0, 0.52, 0, 0.47), ("b", 0.52, 1, 0, 0.47),
            ("c", 0, 0.52, 0.47, 0.935), ("d", 0.52, 1, 0.47, 0.935)]
    },
    "figure_B1": _grid([0, 0.5, 1], [0, 0.5, 1]),
    "figure_C1": {
        "a": _panel((0, 0, 0.49, 0.33)),
        **{letter: _panel((x0, y0, x1, y1), colour_scale=True,
                          footer=(0, 0.953, 1, 1))
           for letter, x0, x1, y0, y1 in [
               ("b", 0.49, 0.92, 0, 0.33),
               ("c", 0, 0.49, 0.33, 0.64),
               ("d", 0.49, 0.92, 0.33, 0.64),
               ("e", 0, 0.49, 0.64, 0.953),
               ("f", 0.49, 0.92, 0.64, 0.953)]},
    },

}


def figure_name(stem):
    name = Path(stem).stem
    match = re.fullmatch(r"figure_(\d+|[ABCS]\d+)", name)
    if not match:
        raise ValueError(f"Unknown figure name: {name}")
    number = match.group(1)
    return "figure_" + (str(int(number)) if number.isdigit() else number)


def png_path(out_stem):
    """Resolve a figure stem or PNG path to its figure output folder."""
    path = Path(out_stem)
    name = figure_name(path)
    folder = path.parent if path.parent.name == name else path.parent / name
    folder.mkdir(parents=True, exist_ok=True)
    return folder / f"{name}.png"


def _trim(image):
    box = ImageChops.difference(image, Image.new("RGB", image.size, "white")).getbbox()
    return image.crop(box) if box else image


def _crop(image, box, trim=True, exclude=()):
    x0, y0, x1, y1 = box
    pixels = (round(x0 * image.width), round(y0 * image.height),
              round(x1 * image.width), round(y1 * image.height))
    part = image.crop(pixels)
    for left, top, right, bottom in exclude:
        region = (round(left * image.width) - pixels[0],
                  round(top * image.height) - pixels[1],
                  round(right * image.width) - pixels[0],
                  round(bottom * image.height) - pixels[1])
        part.paste("white", region)
    return _trim(part) if trim else part


def _stack(parts, gap):
    width = max(part.width for part in parts)
    height = sum(part.height for part in parts) + gap * (len(parts) - 1)
    result = Image.new("RGB", (width, height), "white")
    y = 0
    for part in parts:
        result.paste(part, ((width - part.width) // 2, y))
        y += part.height + gap
    return result


def _with_rank_scale(panel):
    """Repeat C1's exact ten-class colour mapping at a standalone map's height."""
    import matplotlib as mpl
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
    from matplotlib.figure import Figure
    import numpy as np

    # The same colours, boundaries and label as the C1 drawing script.
    cmap = LinearSegmentedColormap.from_list(
        "event_frequency_q1_muted", ["#2ca25f", "#f7fcb9", "#f03b20"])
    bounds = np.round(np.arange(0, 1.0001, 0.1), 1)
    dpi = 500
    with mpl.rc_context({"font.family": "Arial", "font.size": 9}):
        figure = Figure(figsize=(0.88, panel.height / dpi), dpi=dpi)
        FigureCanvasAgg(figure)
        ax = figure.add_axes([0.07, 0.08, 0.20, 0.84])
        bar = ColorbarBase(ax, cmap=cmap,
                          norm=BoundaryNorm(bounds, cmap.N, clip=True),
                          boundaries=bounds, ticks=bounds, alpha=0.78)
        bar.set_label("Relative centroid-concentration rank, R(s)", fontsize=9)
        bar.ax.tick_params(labelsize=9, length=2.5, width=0.6)
        bar.outline.set_linewidth(0.6)
        buffer = io.BytesIO()
        figure.savefig(buffer, format="png", dpi=dpi, facecolor="white")
    buffer.seek(0)
    with Image.open(buffer) as image:
        scale = image.convert("RGB")
    result = Image.new("RGB", (panel.width + scale.width, panel.height), "white")
    result.paste(panel, (0, 0))
    result.paste(scale, (panel.width, 0))
    return result


def export_panels(full_png, *, layout=None, panel_source=None):
    """Export panel regions at native pixel resolution; never alter the full PNG."""
    full_png = Path(full_png)
    name = figure_name(full_png)
    if layout is None:
        layout = PANEL_LAYOUTS[name]
    if not layout:
        return []
    folder = full_png.parent / f"{name}_panels"
    folder.mkdir(parents=True, exist_ok=True)
    written = []
    with Image.open(panel_source or full_png) as source:
        dpi = source.info.get("dpi", (600, 600))
        image = source.convert("RGB")
    if panel_source is not None:
        with Image.open(full_png) as full:
            dpi = full.info.get("dpi", dpi)
    for letter, spec in layout.items():
        keep_width = spec.get("keep_width", False)
        body = _crop(image, spec["box"], trim=not keep_width,
                     exclude=spec.get("exclude", ()))
        if spec.get("colour_scale"):
            body = _with_rank_scale(body)
        parts = []
        if "header" in spec:
            parts.append(_crop(image, spec["header"], trim=not keep_width))
        parts.append(body)
        if "footer" in spec:
            parts.append(_crop(image, spec["footer"]))
        gap = max(16, round(min(body.size) * 0.025))
        panel = _trim(_stack(parts, gap))
        padded = Image.new("RGB", (panel.width + 2 * gap, panel.height + 2 * gap), "white")
        padded.paste(panel, (gap, gap))
        path = folder / f"{name}{letter}.png"
        padded.save(path, dpi=dpi)
        written.append(path)
    return written


def save_figure(fig, out_stem, *, physical_width=None, panel_layout=None, **options):
    """Save the labelled full figure and unlettered standalone panels."""
    path = png_path(out_stem)
    fig.savefig(path, **options)
    if physical_width is not None:
        with Image.open(path) as source:
            image = source.copy()
        dpi = image.width / physical_width
        image.save(path, dpi=(dpi, dpi))
    # A second render supplies the panels. Keep the text artists in place so
    # bbox_inches='tight' produces the same canvas, but omit panel letters.
    from matplotlib.text import Text

    changed = []
    for artist in fig.findobj(match=Text):
        label = artist.get_text()
        match = re.match(r"^\(([a-z])\)(?:\s+|$)", label)
        if match:
            changed.append((artist, label, artist.get_color()))
            if match.end() == len(label):
                artist.set_color("white")
            else:
                artist.set_text(label[match.end():])
    try:
        if changed:
            panel_image = io.BytesIO()
            fig.savefig(panel_image, format="png", **options)
            panel_image.seek(0)
            panels = export_panels(path, layout=panel_layout,
                                   panel_source=panel_image)
        else:
            panels = export_panels(path, layout=panel_layout)
    finally:
        for artist, label, color in changed:
            artist.set_text(label)
            artist.set_color(color)
    return [path, *panels]


def export_powerpoint(source, out_stem):
    """Render a one-slide editable PPTX with desktop PowerPoint on Windows."""
    source = Path(source).resolve(strict=True)
    powershell = shutil.which("powershell.exe") if os.name == "nt" else None
    if powershell is None:
        raise RuntimeError(
            "Figures 1 and 4 require Windows and desktop Microsoft PowerPoint. "
            "Install PowerPoint and the slide fonts (Aptos, Abadi and Cambria Math)."
        )

    # Environment variables carry paths literally, including spaces and quotes.
    # Open the source read-only without a window; close only this presentation.
    script = r'''
$ErrorActionPreference = 'Stop'
[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new($false)
$application = $null
$presentation = $null
$slide = $null
try {
    $application = New-Object -ComObject PowerPoint.Application
    $presentation = $application.Presentations.Open($env:SCORCH_PPTX_SOURCE, -1, 0, 0)
    if ($presentation.Slides.Count -ne 1) {
        throw 'The figure source must contain exactly one slide.'
    }
    $slide = $presentation.Slides.Item(1)
    $slide.Export($env:SCORCH_PPTX_EXPORT, 'PNG', 4500, 2531)
}
catch {
    [Console]::Error.WriteLine($_.Exception.Message)
    exit 1
}
finally {
    if ($null -ne $slide) { [void][Runtime.InteropServices.Marshal]::ReleaseComObject($slide) }
    if ($null -ne $presentation) {
        $presentation.Close()
        [void][Runtime.InteropServices.Marshal]::ReleaseComObject($presentation)
    }
    if ($null -ne $application) { [void][Runtime.InteropServices.Marshal]::ReleaseComObject($application) }
    [GC]::Collect()
    [GC]::WaitForPendingFinalizers()
}
'''
    path = png_path(out_stem)
    temporary = path.resolve().with_name(f".{path.stem}-{uuid4().hex}.png")
    try:
        env = {**os.environ, "SCORCH_PPTX_SOURCE": str(source),
               "SCORCH_PPTX_EXPORT": str(temporary)}
        result = subprocess.run(
            [powershell, "-NoProfile", "-NonInteractive", "-Command", script],
            env=env, capture_output=True, text=True, encoding="utf-8", errors="replace",
            timeout=60, creationflags=subprocess.CREATE_NO_WINDOW,
        )
        if result.returncode != 0:
            detail = result.stderr.strip() or result.stdout.strip()
            raise RuntimeError(f"PowerPoint could not export {source.name}: {detail}")
        with Image.open(temporary) as image:
            if image.size != (4500, 2531):
                raise ValueError(f"Unexpected PowerPoint export size: {image.size}")
            image.verify()
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    return [path, *export_panels(path)]
