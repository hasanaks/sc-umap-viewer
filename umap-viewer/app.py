from shiny.express import ui, render, input
import scanpy as sc
from h5ad import H5AD

ui.input_file("dataset_file", "Select a dataset file (.zip)", accept=[".zip"])


@render.plot(alt="UMAP Plot")
def dataset_plot():
    dataset_file = input.dataset_file()

    if dataset_file:
        dataset = H5AD(dataset_file[0]["datapath"])

        axes = sc.pl.umap(dataset.adata, show=False, title=dataset.name)
        return axes
