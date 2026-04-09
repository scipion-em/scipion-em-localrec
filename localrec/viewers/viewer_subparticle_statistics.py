# ***************************************************************************
# *
# * Authors:   Scipion localrec maintainers
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# ***************************************************************************

import json
import os

from pyworkflow.viewer import DESKTOP_TKINTER, Viewer
from pwem.viewers import EmPlotter

from localrec.protocols.protocol_subparticle_statistics import ProtSubparticleStatistics


class ProtSubparticleStatisticsViewer(Viewer):
    """Visualize histogram outputs from ProtSubparticleStatistics."""

    _label = 'viewer subparticle statistics'
    _targets = [ProtSubparticleStatistics]
    _environments = [DESKTOP_TKINTER]

    def _visualize(self, obj, **args):
        json_path = self.protocol._getExtraPath('histogram.json')
        if not os.path.exists(json_path):
            raise Exception('Histogram output not found: %s' % json_path)

        with open(json_path) as handle:
            data = json.load(handle)

        mode = data.get('mode', '1d')
        if mode == '2d':
            return [self._visualize2D(data)]
        return [self._visualize1D(data)]

    def _visualize1D(self, data):
        bins = data.get('bins', [])
        counts = data.get('counts', [])
        class_ids = data.get('targetClassIds')
        if class_ids is None:
            class_ids = [data.get('targetClass', self.protocol.targetClass.get())]
        class_ids_str = ','.join([str(cid) for cid in class_ids])

        plotter = EmPlotter(windowTitle='Subparticle statistics')
        axis = plotter.createSubPlot('Subparticle occupancy (class ID(s): %s)'
                                     % class_ids_str,
                                     'Subparticles in selected class ID(s)',
                                     'Number of parent particles')
        axis.bar(bins, counts)
        axis.set_xticks(bins)
        return plotter

    def _visualize2D(self, data):
        points = data.get('points', [])
        group1_ids = data.get('targetClassIds', [])
        group2_ids = data.get('targetClass2Ids', [])

        xs = [p['x'] for p in points]
        ys = [p['y'] for p in points]
        freqs = [p['frequency'] for p in points]
        max_freq = max(freqs) if freqs else 1

        sizes = [80.0 + 520.0 * (float(freq) / float(max_freq)) for freq in freqs]

        plotter = EmPlotter(windowTitle='Subparticle statistics')
        axis = plotter.createSubPlot(
            'Subparticle occupancy bubble plot',
            'Subparticles in class ID(s) group 1: %s' % ','.join([str(x) for x in group1_ids]),
            'Subparticles in class ID(s) group 2: %s' % ','.join([str(y) for y in group2_ids]))

        axis.scatter(xs, ys, s=sizes, alpha=0.6)

        for x_value, y_value, freq in zip(xs, ys, freqs):
            axis.text(x_value, y_value, str(freq), ha='center', va='center')

        unique_x = sorted(set(xs))
        unique_y = sorted(set(ys))
        if unique_x:
            axis.set_xticks(unique_x)
        if unique_y:
            axis.set_yticks(unique_y)

        return plotter
