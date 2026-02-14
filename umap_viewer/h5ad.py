"""For handling H5AD data."""

from tempfile import TemporaryDirectory
import scanpy as sc
from pathlib import Path
import zipfile


class H5AD:
    """Holds dataset data."""

    def __init__(self, fpath: str):
        with open(fpath, "rb") as f:
            self.name = Path(f.name).stem
            """Stem of the filepath (e.g., /dir/data123.h5ad will be data123)."""

            self.adata = sc.read_h5ad(f)
            """Holds h5ad data."""


def h5ads_from_zip(zip_path: str) -> list[H5AD]:
    """Load all h5ad files from a ZIP archive.

    :param zip_path: Path to the ZIP archive.
    :return: List of H5AD objects.
    :raises FileNotFoundError: If the file is not a valid ZIP or contains no .h5ad files.
    """

    h5ads = list()

    if not zipfile.is_zipfile(zip_path):
        raise FileNotFoundError(f"The dataset file must be a ZIP file.")

    with TemporaryDirectory() as temp:
        with zipfile.ZipFile(zip_path, "r") as z:
            z.extractall(temp)

        for h5ad_file in Path(temp).glob("*.h5ad"):
            h5ads.append(H5AD(h5ad_file))

    if len(h5ads) == 0:
        raise FileNotFoundError(f"Could not find any .h5ad file in the archive.")

    return h5ads
