[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Leonauna/Educational_CADD/master)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Leonauna/Educational_CADD/blob/master/CADD_pipeline.ipynb)

# Educational CADD Tutorial
An educational tutorial to showcase Computer Aided Drug Discovery (CADD) in the classroom. The interactive Jupyter Notebooks are designed for teaching purposes and lead the user through a simplifed virtual molecular compound screening pipeline.

Our tutorial runs through an example of the protein estrogen receptor (1G50) which has a correlated link to breast cancer. Using a data subset from the ZINC15 database, we apply a random forest classifier to predict the effectiveness of a molecular footprint as a promising drug target.

### [Click here for presentation](https://leonauna.github.io/Educational_CADD/)

# Set up environment

## Setup locally
The easist way to setup an environment for this project is to use [`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) package manager  
`conda env create -f environment.yml`  
To activate the environment `conda activate CADD` 

## Run in binder
You can run the notebook in binder by following this [link](https://mybinder.org/v2/gh/Leonauna/Educational_CADD/master). Note that binder spawn a new container when you launch the notebook so it can take awhile :)

## Run in Colab
You can run the notebook in binder by following this [link](https://colab.research.google.com/github/Leonauna/Educational_CADD/blob/master/CADD_pipeline.ipynb). You may need to install some additional modules when running in colab, you can do so by using '!pip'.

# Contributors
This project was coded during the Life Science Hackathon 2019 in London.\
Contributers are:

[Jonathan Ish-Horowicz](https://github.com/jonathanishhorowicz)\
[Daniel Jiang](https://github.com/WizardOfAus)\
[Léonie Strömich](https://github.com/Leonauna)\
[Tony Yang](https://github.com/tonyyzy) 
