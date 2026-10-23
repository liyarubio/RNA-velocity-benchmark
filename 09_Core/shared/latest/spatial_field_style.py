"""Shared spatial-field display, using the published legacy ED15c parameters."""
import numpy as np
from scipy.spatial import cKDTree

TYPES = ['Inner EPC 1', 'Inner EPC 2', 'Outer EPC 1', 'Outer EPC 2',
         'P-TGC (Sec. P-TGC 1)', 'P-TGC (Sec. P-TGC 2)']
PALETTE = dict(zip(TYPES, ['#A8DADC', '#FFD6A5', '#B7E4C7', '#FFADAD', '#CDB4DB', '#DDBEA9']))
PARAMETERS = dict(aspect='equal', grid=32, neighbors=30, smoothing_factor=.75,
                  min_mass=.2, padding_fraction=.05, density=.9, linewidth=.43,
                  arrowsize=.53, arrow_color='#292929', point_size=2.8, alpha=.9,
                  background_color='#d2d2d2')

def spatial_field(ax, xy, labels, vectors=None, size=None, full_palette=None):
    """Draw saved coordinates/vectors without fitting or projecting velocities."""
    xy, labels = np.asarray(xy), np.asarray(labels)
    p = PARAMETERS
    target = np.isin(labels, TYPES)
    groups = ([(labels == k, color) for k, color in full_palette.items()]
              if full_palette is not None else
              [(~target, p['background_color'])] + [(labels == k, color) for k, color in PALETTE.items()])
    for keep, color in groups:
        ax.scatter(xy[keep, 0], xy[keep, 1], s=p['point_size'] if size is None else size,
                   c=color, edgecolors='none', alpha=p['alpha'], rasterized=False, zorder=1)
    span = np.ptp(xy, axis=0)
    lo, hi = xy.min(0)-span*p['padding_fraction'], xy.max(0)+span*p['padding_fraction']
    if vectors is not None:
        v = np.asarray(vectors)
        ok = np.isfinite(v).all(1) & (np.linalg.norm(v, axis=1) > 0)
        assert ok.sum() >= 2, 'Insufficient saved valid directions'
        grid = [np.linspace(lo[i], hi[i], p['grid']) for i in range(2)]
        xx, yy = np.meshgrid(*grid)
        pts = np.c_[xx.ravel(), yy.ravel()]
        ds, idx = cKDTree(xy[ok]).query(pts, k=min(p['neighbors'], ok.sum()))
        scale = np.mean(span)/p['grid']*p['smoothing_factor']
        weights = np.exp(-.5*(ds/max(scale, 1e-8))**2)
        mass = weights.sum(1)
        field = (v[ok][idx]*weights[:, :, None]).sum(1)/np.maximum(mass[:, None], 1e-10)
        field[mass < p['min_mass']] = np.nan
        ax.streamplot(*grid, np.ma.masked_invalid(field[:, 0].reshape(xx.shape)),
                      np.ma.masked_invalid(field[:, 1].reshape(xx.shape)),
                      density=p['density'], color=p['arrow_color'], linewidth=p['linewidth'],
                      arrowsize=p['arrowsize'], zorder=2)
    ax.set(xlim=(lo[0], hi[0]), ylim=(lo[1], hi[1]))
    ax.set_aspect(p['aspect'])
    ax.axis('off')
