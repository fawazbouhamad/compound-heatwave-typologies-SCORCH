# SPDX-License-Identifier: GPL-3.0-only
"""SCORCH: Spatiotemporal Classification of Regional Compound Heatwaves.

The package holds the calculations of the five-stage SCORCH workflow:

1. ``scorch.heatwaves``   local heatwaves at each grid cell
2. ``scorch.regional``    regionally extensive heatwave days
3. ``scorch.clustering``  DBSCAN clusters of heatwave cells on each selected day
4. ``scorch.ellipses``    PCA ellipses and temperature-weighted centroids
5. ``scorch.events``      compound events and their four types (Table 1)

``scorch.stages`` runs the stages on a prepared NetCDF file, ``scorch.config``
holds the study settings, ``scorch.io`` reads and writes files, and
``scorch.demo`` builds the small synthetic example. ``scorch.cli`` is the
``scorch`` command.
"""
from __future__ import annotations

__version__ = "1.0.0"
__all__ = ["__version__"]
