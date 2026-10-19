import unittest
import numpy as np
from scipy import sparse
from evaluate_reference import aligned_scores


class ReferenceTests(unittest.TestCase):
    def test_reorders_cells_and_keeps_missing(self):
        result = aligned_scores(['b', 'a'], [[1, 0], [0, 0]], [[0, 1], [1, 0]],
                                ['a', 'b', 'c'], [[0, 0], [1, 0], [2, 0]], [[1, 0], [0, -1], [1, 1]])
        np.testing.assert_allclose(result.cosine[:2], [1, -1])
        self.assertEqual(result.status[2], 'missing_prediction')
        self.assertTrue(np.isnan(result.cosine[2]))

    def test_invalid_vectors_and_coordinates(self):
        result = aligned_scores(['a', 'b', 'c'], [[0, 0], [0, 0], [1, 1]],
                                [[0, 0], [np.nan, 1], [1, 0]], ['a', 'b', 'c'],
                                [[0, 0], [0, 0], [0, 0]], [[1, 0]] * 3)
        self.assertEqual(result.status.tolist(), ['zero_vector', 'nonfinite_velocity', 'coordinate_mismatch'])
        self.assertTrue(result.cosine.isna().all())

    def test_sparse_vectors(self):
        value = sparse.csr_matrix([[1., 0.]])
        result = aligned_scores(['a'], [[0, 0]], value, ['a'], [[0, 0]], value)
        self.assertEqual(result.cosine[0], 1)

    def test_duplicate_ids_rejected(self):
        with self.assertRaisesRegex(ValueError, 'unique'):
            aligned_scores(['a', 'a'], [[0, 0]] * 2, [[1, 1]] * 2,
                           ['a'], [[0, 0]], [[1, 1]])


if __name__ == '__main__':
    unittest.main()
