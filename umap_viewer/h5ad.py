from typing import IO
import scanpy as sc
from pathlib import Path
import zipfile


class H5AD:
    # reads an .h5ad file from a file-like object
    def __init__(self, file: IO[bytes]):
        self.name = str(Path(file.name).with_suffix(""))
        self.adata = sc.read_h5ad(file)


# reads and returns a list[H5AD] from a zip archive
def h5ads_from_zip(zip_path: str) -> list[H5AD]:
    h5ads = list()

    with zipfile.ZipFile(zip_path, "r") as z:
        for file in z.namelist():
            if file.endswith(".h5ad"):
                with z.open(file) as dataset:
                    h5ads.append(H5AD(dataset))

    if len(h5ads) == 0:
        raise FileNotFoundError(
            f"Could not find any .h5ad file in the archive {zip_path}."
        )

    return h5ads
