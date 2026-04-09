# **************************************************************************
# *
# * Authors:   Scipion localrec maintainers
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# **************************************************************************

import csv
import json
import os

from pyworkflow import VERSION_3_0
from pyworkflow.protocol.params import PointerParam, IntParam, StringParam
from pwem.protocols import ProtParticles
from localrec.utils import get_subparticle_provenance


class ProtSubparticleStatistics(ProtParticles):
    """Compute class-occupancy statistics for subparticles.

    1D mode (group 1 only):
        - X axis categories: number of subparticles assigned to class ID(s)
          per original particle (0..N).
        - Y axis: number of original particles in each category.

    2D mode (group 1 + group 2):
        - X axis categories: number of subparticles assigned to class ID(s)
          in group 1 per original particle.
        - Y axis categories: number of subparticles assigned to class ID(s)
          in group 2 per original particle.
        - Frequency: number of original particles with each (X, Y) pair.
    """

    _label = 'subparticle statistics'
    _lastUpdateVersion = VERSION_3_0

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputClassification', PointerParam,
                      pointerClass='SetOfClasses3D, SetOfClasses2D',
                      important=True,
                      label='Input classification',
                      help='Select a classification output where each class '
                           'contains subparticles.')

        form.addParam('targetClass', StringParam, default='1',
                      label='Class ID(s)',
                      help='Class ID or comma-separated class IDs used to '
                           'compute occupancy histogram. Example: 1,2,5')

        form.addParam('targetClass2', StringParam, default='',
                      label='Class ID(s) group 2 (optional)',
                      help='Optional second class group (single ID or '
                           'comma-separated IDs). If provided, a 2D histogram '
                           'is computed and exported.')

        form.addParam('maxCount', IntParam, default=-1,
                      label='Maximum count N (optional)',
                      help='If > 0, force histogram limits up to N. If <= 0, '
                           'limits are automatically inferred from '
                           'classification data.')

        form.addParallelSection(threads=0, mpi=0)

    def _insertAllSteps(self):
        self._insertFunctionStep('computeStatisticsStep')

    def computeStatisticsStep(self):
        classification = self.inputClassification.get()
        group1_ids = self._parseClassIds(self.targetClass.get(),
                                         'Class ID(s)',
                                         allow_empty=False)
        group2_ids = self._parseClassIds(self.targetClass2.get(),
                                         'Class ID(s) group 2 (optional)',
                                         allow_empty=True)

        per_particle_counts = {}
        provenance_summary = {'withProvenance': 0, 'missingProvenance': 0}

        for parent_id, class_id, has_provenance in self._iterAssignments(classification):
            if parent_id not in per_particle_counts:
                per_particle_counts[parent_id] = [0, 0]
            if has_provenance:
                provenance_summary['withProvenance'] += 1
            else:
                provenance_summary['missingProvenance'] += 1

            if class_id in group1_ids:
                per_particle_counts[parent_id][0] += 1

            if group2_ids and class_id in group2_ids:
                per_particle_counts[parent_id][1] += 1

        if not per_particle_counts:
            raise Exception('No subparticle assignments found in the selected '
                            'classification input.')
        if provenance_summary['missingProvenance'] > 0:
            self.warning('Subparticle provenance fields '
                         '(_symmetryGroup/_symmetryOperatorId) were missing in '
                         '%d items. Continuing statistics computation.'
                         % provenance_summary['missingProvenance'])

        if group2_ids:
            self._write2DHistogram(per_particle_counts, group1_ids, group2_ids,
                                   provenance_summary)
        else:
            self._write1DHistogram(per_particle_counts, group1_ids,
                                   provenance_summary)

    def _write1DHistogram(self, per_particle_counts, group1_ids, provenance_summary):
        counts_group1 = [count_pair[0] for count_pair in per_particle_counts.values()]

        auto_n = max(counts_group1)
        n = self._resolveMaxCount(auto_n)

        histogram = {k: 0 for k in range(n + 1)}
        overflow = 0
        for count in counts_group1:
            if count <= n:
                histogram[count] += 1
            else:
                overflow += 1

        if overflow:
            self.warning('%d particles exceeded maximum bin N=%d and were '
                         'grouped in overflow.' % (overflow, n))

        bins = list(range(n + 1))
        data = {
            'mode': '1d',
            'targetClass': self.targetClass.get(),
            'targetClassIds': sorted(group1_ids),
            'maxCount': n,
            'overflow': overflow,
            'totalParticles': len(per_particle_counts),
            'provenance': provenance_summary,
            'bins': bins,
            'counts': [histogram[k] for k in bins]
        }

        json_path, csv_path = self._getHistogramPaths()
        with open(json_path, 'w') as handle:
            json.dump(data, handle, indent=2)

        with open(csv_path, 'w') as handle:
            writer = csv.writer(handle)
            writer.writerow(['category_k', 'particle_count',
                             'with_provenance', 'missing_provenance'])
            for k in bins:
                writer.writerow([k, histogram[k],
                                 provenance_summary['withProvenance'],
                                 provenance_summary['missingProvenance']])

    def _write2DHistogram(self, per_particle_counts, group1_ids, group2_ids,
                          provenance_summary):
        pairs = [tuple(count_pair) for count_pair in per_particle_counts.values()]
        x_auto = max(pair[0] for pair in pairs)
        y_auto = max(pair[1] for pair in pairs)
        n = self._resolveMaxCount(max(x_auto, y_auto))

        matrix = {}
        overflow = 0
        for x_count, y_count in pairs:
            if x_count <= n and y_count <= n:
                key = (x_count, y_count)
                matrix[key] = matrix.get(key, 0) + 1
            else:
                overflow += 1

        points = []
        for key in sorted(matrix.keys()):
            points.append({'x': key[0], 'y': key[1], 'frequency': matrix[key]})

        if overflow:
            self.warning('%d particles exceeded maximum limits N=%d and were '
                         'grouped in overflow.' % (overflow, n))

        data = {
            'mode': '2d',
            'targetClass': self.targetClass.get(),
            'targetClassIds': sorted(group1_ids),
            'targetClass2': self.targetClass2.get(),
            'targetClass2Ids': sorted(group2_ids),
            'maxCount': n,
            'overflow': overflow,
            'totalParticles': len(per_particle_counts),
            'provenance': provenance_summary,
            'points': points
        }

        json_path, csv_path = self._getHistogramPaths()
        with open(json_path, 'w') as handle:
            json.dump(data, handle, indent=2)

        with open(csv_path, 'w') as handle:
            writer = csv.writer(handle)
            writer.writerow(['x_count_group1', 'y_count_group2',
                             'particle_frequency', 'with_provenance',
                             'missing_provenance'])
            for point in points:
                writer.writerow([point['x'], point['y'], point['frequency'],
                                 provenance_summary['withProvenance'],
                                 provenance_summary['missingProvenance']])

    def _resolveMaxCount(self, auto_n):
        user_n = int(self.maxCount.get())
        n = user_n if user_n > 0 else auto_n
        if n < auto_n:
            self.warning('Provided maxCount (%d) is smaller than observed '
                         'maximum (%d). Using observed maximum.' % (n, auto_n))
            n = auto_n
        return n

    def _getHistogramPaths(self):
        return self._getExtraPath('histogram.json'), self._getExtraPath('histogram.csv')

    def _iterAssignments(self, classification):
        """Yield tuples: (parent_particle_id, class_id)."""
        for cls in classification:
            class_id = self._getClassId(cls)
            if hasattr(cls, 'iterItems'):
                cls_items = cls.iterItems()
            else:
                cls_items = cls
            for item in cls_items:
                parent_id = self._getParentParticleId(item)
                sym_group, operator_id = get_subparticle_provenance(item)
                has_provenance = sym_group is not None and operator_id is not None
                yield parent_id, class_id, has_provenance

    @staticmethod
    def _getClassId(cls):
        for attr in ('_objId', 'getObjId', 'id'):
            if hasattr(cls, attr):
                value = getattr(cls, attr)
                value = value() if callable(value) else value
                if value is not None:
                    return int(value)
        raise Exception('Could not read class id from classification class.')

    @staticmethod
    def _getParentParticleId(subparticle):
        if hasattr(subparticle, 'getCoordinate'):
            coord = subparticle.getCoordinate()
            if coord is not None and hasattr(coord, '_micId'):
                return int(coord._micId)

        if hasattr(subparticle, '_micId'):
            return int(subparticle._micId)

        raise Exception('Could not read parent particle id for subparticle.')

    def _parseClassIds(self, raw_value, param_label, allow_empty):
        if raw_value is None:
            raw_value = ''
        raw_items = [x.strip() for x in str(raw_value).split(',')]
        class_ids = set()

        for item in raw_items:
            if item == '':
                continue
            try:
                value = int(item)
            except Exception:
                raise Exception('%s must be integer values separated by '
                                'commas. Example: 1,2,5' % param_label)
            if value < 1:
                raise Exception('%s must be >= 1.' % param_label)
            class_ids.add(value)

        if not class_ids and not allow_empty:
            raise Exception('%s can not be empty.' % param_label)

        return class_ids

    def _validate(self):
        errors = []
        try:
            group1 = self._parseClassIds(self.targetClass.get(),
                                         'Class ID(s)',
                                         allow_empty=False)
            group2 = self._parseClassIds(self.targetClass2.get(),
                                         'Class ID(s) group 2 (optional)',
                                         allow_empty=True)
            overlap = group1.intersection(group2)
            if overlap:
                errors.append('Class ID(s) and Class ID(s) group 2 (optional) '
                              'must not overlap. Overlap: %s' % sorted(overlap))
        except Exception as ex:
            errors.append(str(ex))
        return errors

    def _summary(self):
        summary = []
        json_path, csv_path = self._getHistogramPaths()
        if os.path.exists(json_path):
            group2 = self._parseClassIds(self.targetClass2.get(),
                                         'Class ID(s) group 2 (optional)',
                                         allow_empty=True)
            if group2:
                summary.append('2D occupancy histogram computed for class ID '
                               'groups: [%s] vs [%s].'
                               % (self.targetClass.get(), self.targetClass2.get()))
            else:
                summary.append('Occupancy histogram computed for class ID(s): %s.'
                               % self.targetClass.get())
            summary.append('Results: %s' % json_path)
            summary.append('Results: %s' % csv_path)
        return summary

    def _methods(self):
        return []
