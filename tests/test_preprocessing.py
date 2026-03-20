from __future__ import annotations

import unittest

import pandas as pd

from mofapipeline.preprocessing import filter_matrix_by_missingness, impute_matrix


class PreprocessingTestCase(unittest.TestCase):
    def test_filter_matrix_by_missingness_removes_sparse_rows_and_columns(self) -> None:
        matrix = pd.DataFrame(
            {
                "s1": [1.0, None, 2.0],
                "s2": [2.0, None, None],
                "s3": [3.0, 1.0, None],
            },
            index=["f1", "f2", "f3"],
        )
        filtered = filter_matrix_by_missingness(matrix, feature_threshold=0.5, sample_threshold=0.5)
        self.assertListEqual(filtered.index.tolist(), ["f1"])
        self.assertListEqual(filtered.columns.tolist(), ["s1", "s3"])

    def test_impute_matrix_uses_feature_median(self) -> None:
        matrix = pd.DataFrame(
            {"s1": [1.0, 10.0], "s2": [None, 30.0], "s3": [3.0, None]},
            index=["f1", "f2"],
        )
        imputed = impute_matrix(matrix)
        self.assertEqual(imputed.loc["f1", "s2"], 2.0)
        self.assertEqual(imputed.loc["f2", "s3"], 20.0)


if __name__ == "__main__":
    unittest.main()
