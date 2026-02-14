from tempfile import TemporaryDirectory
import scanpy as sc
from pathlib import Path
import zipfile


class H5AD:
    # reads an .h5ad file into AnnData
    def __init__(self, fpath: str):
        with open(fpath, "rb") as f:
            self.name = Path(f.name).stem
            self.adata = sc.read_h5ad(f)


# reads and returns a list[H5AD] from a zip archive
def h5ads_from_zip(zip_path: str) -> list[H5AD]:
    h5ads = list()

    with TemporaryDirectory() as temp:
        with zipfile.ZipFile(zip_path, "r") as z:
            z.extractall(temp)

        for h5ad_file in Path(temp).glob("*.h5ad"):
            h5ads.append(H5AD(h5ad_file))

    if len(h5ads) == 0:
        raise FileNotFoundError(
            f"Could not find any .h5ad file in the archive {zip_path}."
        )

    return h5ads
