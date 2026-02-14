from shiny.express import ui, render, input
from shiny import reactive, req
import scanpy as sc
import h5ad

ui.page_opts(title="UMAP Viewer", fillable=True)

with ui.sidebar():
    ui.input_file("dataset_archive", "Select a dataset file (.zip)", accept=[".zip"])
    ui.input_select("selected_dataset", "Select a dataset", [])


with ui.card():

    @render.plot(alt="UMAP Plot")
    def dataset_plot():
        dataset_name = req(input.selected_dataset())
        dataset_list = req(datasets())

        for dataset in dataset_list:
            if dataset.name == dataset_name:
                data = dataset

        axes = sc.pl.umap(data.adata, show=False, title=data.name)
        return axes


@reactive.calc
def datasets() -> list[h5ad.H5AD]:
    archive = req(input.dataset_archive())
    return h5ad.h5ads_from_zip(archive[0]["datapath"])


@reactive.effect
def _():
    ui.update_select(
        "selected_dataset",
        choices=[dataset.name for dataset in datasets()],
    )
