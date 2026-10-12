"""Shared layout functions for the three representative-event figures."""
import matplotlib.pyplot as plt
import figure_style as Q
from figure_exports import save_figure
import snapshot_maps as SN

DPI = 600
TIME_FS = 15
def check_event(compound_events, dates, expected_event, expected_type):
    """Stop unless all dates belong to the expected event and type."""
    for date in dates:
        row = compound_events[(compound_events["start_date"] <= date)
                              & (compound_events["end_date"] >= date)]
        assert len(row) == 1, date
        row = row.iloc[0]
        assert int(row["compound_event_id"]) == expected_event, (date, row)
        assert int(row["event_type"]) == expected_type, (date, row)

def _render(ax, cluster_cells, ellipses, date):
    day_cells = cluster_cells[cluster_cells["date"] == date]
    day_ellipses = ellipses[ellipses["date"] == date]
    SN.render_day(ax, day_cells, day_ellipses)
    Q.date_label(ax, date)

def _grid(nrows, ncols, figsize, top, bottom, hspace=0.10):
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False,
                             subplot_kw={"projection": SN.PROJ})
    fig.subplots_adjust(top=top, bottom=bottom, left=0.055, right=0.985,
                        wspace=0.05, hspace=hspace)
    return fig, axes

def _time_label(fig, ax, text):
    b = ax.get_position()
    fig.text((b.x0 + b.x1) / 2.0, b.y1 + 0.012,
             rf"$\mathbf{{{text}}}$", ha="center", va="bottom",
             fontsize=TIME_FS, color="black")

def _save(fig, out_stem):
    written = save_figure(fig, out_stem, dpi=DPI)
    plt.close(fig)
    return written
