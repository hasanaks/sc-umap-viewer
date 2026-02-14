from umap_viewer import h5ad
from zipfile import BadZipFile
import os
import pytest


def test_valid_zip_valid_h5ad():
    path = os.path.join("data", "pbmc3k.zip")
    datasets = h5ad.h5ads_from_zip(path)

    assert len(datasets) == 1
    assert datasets[0].adata is not None
    assert datasets[0].name == "pbmc3k"


def test_valid_zip_no_h5ad():
    with pytest.raises(FileNotFoundError):
        path = os.path.join("tests", "empty.zip")
        h5ad.h5ads_from_zip(path)


def test_non_zip():
    with pytest.raises(BadZipFile):
        path = os.path.join("tests", "empty.txt")
        h5ad.h5ads_from_zip(path)
