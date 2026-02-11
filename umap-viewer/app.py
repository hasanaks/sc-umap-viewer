from shiny.express import ui, render, input
from zipfile import ZipFile

ui.input_file("dataset", "Select a dataset file (.zip)", accept=[".zip"])

# @render.plot()
# def dataset_plot():
#     pass

@render.text
def file_contents():
    dataset = input.dataset()

    if dataset:
        print("opened:", dataset)

        with ZipFile(dataset[0]["datapath"]) as dataset_zip:
            for fpath in dataset_zip.namelist():
                if fpath.endswith(".h5ad"):
                    with dataset_zip.open(fpath) as dataset_file:
                        pass # todo: generate umap plot
