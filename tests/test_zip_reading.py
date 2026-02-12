import unittest
from umap_viewer.h5ad import H5AD
import zipfile
import os


class TestZipReading(unittest.TestCase):

    def test_valid_zip_valid_h5ad(self):
        path = os.path.join("data", "pbmc3k.zip")
        dataset = H5AD(path)
        self.assertIsNotNone(dataset.adata)
        self.assertEqual(dataset.name, "pbmc3k")

    def test_valid_zip_no_h5ad(self):
        with self.assertRaises(FileNotFoundError):
            path = os.path.join("tests", "empty.zip")
            H5AD(path)

    def test_non_zip(self):
        with self.assertRaises(zipfile.BadZipFile):
            path = os.path.join("tests", "empty.txt")
            H5AD(path)


if __name__ == "__main__":
    unittest.main()
