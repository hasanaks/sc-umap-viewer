from shiny.express import ui, render, input
import zipfile
import scanpy as sc
from h5ad import H5AD

ui.input_file("dataset", "Select a dataset file (.zip)", accept=[".zip"])


@render.plot(alt="UMAP Plot")
def dataset_plot():
    dataset = input.dataset()

    if dataset:
        dataset = H5AD(dataset[0]["datapath"])

        axes = sc.pl.umap(dataset.adata, show=False)
        return axes
