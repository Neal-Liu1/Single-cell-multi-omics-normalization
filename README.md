Contains the code for testing fast-RUVIII on a few CITE-seq datasets, as well as the reproduction of [Ahlmann-Eltze & Huber's Nature Methods benchmarking study](https://www.nature.com/articles/s41592-023-01814-1). 
A package for fastRUVIII (still in development) is available [here](https://github.com/Neal-Liu1/fastRUVIII)

## Main analysis files for reproducing the figures:
- **Paper panels**: R markdown of the single batch normalization benchmark of [Stuart et al's 2019 bone marrow CITE-seq](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419305598%3Fshowall%3Dtrue)
- **Paper panels 2**: R markdown of the integration benchmark of Stuart et al & [Triana et al's 2021 bone marrow CITE-seq](https://www.nature.com/articles/s41590-021-01059-0)
- **Paper panels 3**: R markdown for benchmark of the [Kotliarov et al's CITE-seq single cell data](https://www.nature.com/articles/s41591-020-0769-8) of baseline PBMC samples from 20 healthy individuals (10 high and 10 low responders)
vaccinated with influenza pandemic H1N1 and seasonal vaccines in 2009.
- **totalVI_runs.ipynb**: The python code (as a jupyter notebook) for running totalVI on all the datasets.
- **Sc_helper_functions.R**: R file of all the helper functions.
- **Ahlmann-Eltze & Huber 2023 consistency.Rmd**: R markdown that acts as a scaffold for the consistency benchmarking. Sets up file paths, params, downloads data, invokes SLURM jobs in the slurm_benchmarking_scripts folder in the right sequence, and collects final results into a final tsv.
- **Ahlmann-Eltze & Huber 2023 downsampling.Rmd**: R markdown that acts as a scaffold for the downsampling benchmarking.
- **Ahlmann-Eltze & Huber 2023 simulation.Rmd**: R markdown that acts as a scaffold for the simulation benchmarking.
- **Ahlmann-Eltze & Huber 2023 plots.Rmd**: All the associated plots in one R markdown.
- **slurm_benchmarking_scripts**: The actual benchmarking code run on the HPC. Contains both shell and R scripts adapted from Ahlmann-Eltze's paper to run on WEHI's SLURM array. 

<img src="./plots/transformgampoi headline.jpg">
