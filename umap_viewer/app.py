from shiny.express import ui, render, input
from shiny import reactive, req
import scanpy as sc
from h5ad import H5AD

ui.page_opts(title="UMAP Viewer", fillable=True)

@reactive.calc
def dataset():
    file = req(input.dataset_file())
    return H5AD(file[0]["datapath"])

with ui.sidebar():
    ui.input_file("dataset_file", "Select a dataset file (.zip)", accept=[".zip"])


with ui.card():
    @render.plot(alt="UMAP Plot")
    def dataset_plot():
        if data := dataset():
            axes = sc.pl.umap(data.adata, show=False, title=data.name)
            return axes
