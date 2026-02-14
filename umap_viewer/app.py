from shiny.express import ui, render, input, expressify
from shiny import reactive, req
import scanpy as sc
import h5ad

ui.page_opts(title="UMAP Viewer", fillable=True)
datasets = reactive.value(None)

with ui.sidebar():
    ui.input_file("dataset_archive", "Select a dataset file (.zip)", accept=[".zip"])
    ui.input_select("selected_dataset", "Select a dataset", [])


with ui.card():

    @render.plot(alt="UMAP Plot")
    def dataset_plot():
        dataset_name = req(input.selected_dataset())
        dataset_list = req(datasets.get())

        for dataset in dataset_list:
            if dataset.name == dataset_name:
                data = dataset

        axes = sc.pl.umap(data.adata, show=False, title=data.name)
        return axes


@reactive.effect
def read_datasets() -> list[h5ad.H5AD]:
    archive = req(input.dataset_archive())

    try:
        datasets.set(h5ad.h5ads_from_zip(archive[0]["datapath"]))
    except Exception as e:
        m = ui.modal(
            str(e)
        )
        ui.modal_show(m)
        
        datasets.set(None)


@reactive.effect
def _():
    dataset_list = datasets.get()

    if dataset_list:
        choices = [""] + [dataset.name for dataset in dataset_list]
    else:
        choices = [""]

    ui.update_select(
        "selected_dataset",
        choices=choices,
    )
