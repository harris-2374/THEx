# Tree House Explorer

Tree House Explorer (THEx) is a novel phylogenomic genome browser that hosts a growing collection of interactive dashboards. It was developed using [Plotly's](https://plotly.com/) open-source framework [Dash](https://dash.plotly.com/) for building extremely powerful and customizable analytical web applications.

## Installation

1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
   - Windows 10 installation [documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html)
   - Mac OS installation [documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html)
   - Linux installation [documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
2. Add required channels
    > $ conda config --add channels bioconda  
    > $ conda config --add channels conda-forge
3. Create a new conda environment
    > $ conda create -n thex_env
4. Activate conda environment
    > $ conda activate thex_env
5. Install THEx
    > $ conda install thex -c ajharris_2374

## Documentation

More information and details for THEx and THExBuilder can be found in the [Wiki](https://github.com/harris-2374/THEx/wiki) tab.

## Citation

If you use THEx in your work, please cite:

  > Andrew J Harris, Nicole M Foley, Tiffani L Williams, William J Murphy, Tree House Explorer: A Novel Genome Browser for Phylogenomics, _Molecular Biology and Evolution_, Volume 39, Issue 6, June 2022, msac130, <https://doi.org/10.1093/molbev/msac130>
