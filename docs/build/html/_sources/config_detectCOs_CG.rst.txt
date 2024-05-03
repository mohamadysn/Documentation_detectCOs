config_detectCOs_CG.yaml
=========================

This is a `config_detectCOs_CG.yaml` file used for creating a Conda environment named `detectCOs`. It specifies the channels where Conda will search for packages (`conda-forge`, `bioconda`, `defaults`, `atgc-montpellier`, and `r`) and lists the dependencies required for the environment.

**Code :**

.. code-block:: yaml

    name: detectCOs
    channels:
    - conda-forge
    - bioconda
    - defaults
    - atgc-montpellier
    - r
    dependencies:
    - python=3.9
    - pip=23.2
    - pip:
        - hues==0.2.2
        - scipy==1.11.1
        - pyyaml==6.0.1
    - r-base=4.2 #changed from 4.3.1 to resolve conflict
    - r-dplyr=1.1.3
    - r-stringr=1.5.0  
    - r-yaml=2.3.7
    - r-data.table=1.14.8
    - r-unix=1.5.4

    
..


This file defines the environment requirements for running the `detectCOs` tool, ensuring that all necessary dependencies are installed in a consistent and reproducible manner.

