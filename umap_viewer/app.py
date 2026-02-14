from shiny.express import ui, render, input
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
        """Plots the selected dataset."""

        dataset_name = req(input.selected_dataset())
        dataset_list = req(datasets.get())

        for dataset in dataset_list:
            if dataset.name == dataset_name:
                data = dataset

        axes = sc.pl.umap(data.adata, show=False, title=data.name)
        return axes


@reactive.effect
def read_datasets():
    """Reads datasets from the uplodad ZIP archive.

    Loads all .h5ad files from the user-selected ZIP archive and stores them
    in the reactive datasets value. Displays an error modal if loading fails."""

    archive = req(input.dataset_archive())

    try:
        datasets.set(h5ad.h5ads_from_zip(archive[0]["datapath"]))
    except Exception as e:
        m = ui.modal(str(e))
        ui.modal_show(m)

        datasets.set(None)


@reactive.effect
def dataset_choices():
    """Updates the list of choices of the selection UI."""

    dataset_list = datasets.get()

    if dataset_list:
        choices = [""] + [dataset.name for dataset in dataset_list]
    else:
        choices = [""]

    ui.update_select(
        "selected_dataset",
        choices=choices,
    )
