import zipfile
import scanpy as sc


class H5AD:
    @staticmethod
    def _find_dataset(archive_path):
        with zipfile.ZipFile(archive_path, "r") as dataset_zip:
            for path in dataset_zip.namelist():
                if path.endswith(".h5ad"):
                    return dataset_zip.open(path, "r")

    # reads an .h5ad file from a zip archive
    def __init__(self, archive_path):
        with H5AD._find_dataset(archive_path) as dataset_file:
            self.adata = sc.read_h5ad(dataset_file)
