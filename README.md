# Roller Coaster Dynamics
This repository contains the interactive reproduction of [Kaschke and Cartarius (2026)](https://doi.org/10.1119/5.0237959) (their arXiv submission can be found [here](https://arxiv.org/abs/2603.27230)) using [Pluto notebook](https://plutojl.org/) and [Julia](https://julialang.org/). To run this notebook locally, please follow the instructions below. In case of errors or doubts, feel free to contact me.

> [!IMPORTANT]
> This repository defines a custom package `PlutoCSSMacros` that I created to make the presentation look good, and it is required to run this Pluto Notebook. Omitting it will result in errors. So clone the entire repo.

## Instructions to run the notebook
Following are the instructions to run the notebook in your local machine:
1. Clone this repository to your local machine.
2. Install Julia from `julialang.org` and Open the terminal.
3. `cd \path\to\cloned\dir`
4. Run `julia --project=.` (This will open the Julia REPL which will look similar to the python REPL if you are familar with that.)
5. Type `import Pluto; Pluto.run()` (This will open a local web-app at `http://localhost:1234/` which might seem familiar if you have used Jupyter notebook before.)
6. Open the `roller_coaster.jl` notebook, *et voilà!*

**Note:** During your first run, you will need to instantiate the packages before you import and run Pluto, so that your environment can fetch the required packages from the web. To do that, press `]` while inside the Julia REPL. This will open the package mode, make sure it says `(JournalClub) pkg>`. Then just type `instantiate` and press Enter to let Julia download the requirements.

You only need to do this the first time between steps 4 and 5, thereafter you can proceed to step 5 directly.

Additionally, when you first run the notebook, it takes a couple minutes to execute all the cells, the status of which can be noted in the `Status` section at the bottom right. There are a couple cells (particularly the ones with interactive elements) that may display error on the first run; do not fret, just execute the cells with the switches again (which should be right above the cells with errors).

> [!WARNING]
> The colours of the notebook are optimised for dark mode, so they may not work properly if your system default is light mode. To change the colours, find the cell with `PlutoCSSMacros.set_colors()` under the *Presentation Helpers* section, and comment it out or change the colors as per your wish (its defintion can be seen upon selection in the `Live docs` at the bottom right).