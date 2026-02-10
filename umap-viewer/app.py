from shiny.express import ui, render, input

ui.input_file("dataset", "Select a dataset file (.zip)", accept=[".zip"])

@render.text
def txt():
    return input.dataset()
