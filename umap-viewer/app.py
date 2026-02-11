from shiny.express import ui, render, input
import zipfile
import scanpy as sc

ui.input_file("dataset", "Select a dataset file (.zip)", accept=[".zip"])


# reads and returns an anndata obj of .h5ad file from a .zip file
def read_dataset(archive_path):
    with zipfile.ZipFile(archive_path, "r") as dataset_zip:
        for path in dataset_zip.namelist():
            if path.endswith(".h5ad"):
                with dataset_zip.open(path, "r") as file:
                    return sc.read_h5ad(file)

@render.plot()
def dataset_plot():
    dataset = input.dataset()

    if dataset:
        print("opened:", dataset)
        adata = read_dataset(dataset[0]["datapath"])
        print(adata)