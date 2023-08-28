# Genetic identification of the common bristlenose

##### Rupert A. Collins // August 2023

![_Ancistrus_ sp.(3)](assets/Ancistrus_sp(3).jpg)


### BACKGROUND

This README file describes the software requirements to run the `ancistr` analysis on an Ubuntu system. The R script `scripts/ancistr.R` searches the [NCBI GenBank](https://www.ncbi.nlm.nih.gov/nuccore) nucleotide database for _Ancistrus_ sequence data (COI gene), and then runs a phylogenetic and sequence similarity analysis to estimate the closest relatives of the common bristlenose catfish [_Ancistrus_ sp. (3)](https://www.planetcatfish.com/ancistrus_cf_cirrhosus). The [Quarto](https://quarto.org/) script `ancistr.qmd` executes this R script and generates an html output report `ancistr.html` (not committed to git).


### INSTALLATION

Install prerequisites from Ubuntu repositories, and install Quarto:
```bash
sudo apt update
sudo apt install make build-essential libssl-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncursesw5-dev xz-utils tk-dev libxml2-dev libxmlsec1-dev libffi-dev liblzma-dev git mafft
# see more recent versions of quarto at https://github.com/quarto-dev/quarto-cli/releases/latest or https://quarto.org/docs/get-started/
wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.3.450/quarto-1.3.450-linux-amd64.deb
sudo apt install ./quarto-1.3.450-linux-amd64.deb
```

Install R (may take up to ten minutes):

```bash
# clone R installer
git clone https://github.com/jcrodriguez1989/renv-installer.git ~/.renv
# define environment variable RENV_ROOT
echo 'export RENV_ROOT="$HOME/.renv"' >> ~/.bashrc
echo 'export PATH="$RENV_ROOT/bin:$PATH"' >> ~/.bashrc
# add renv init to your shell
echo -e 'if command -v renv 1>/dev/null 2>&1; then\n  eval "$(renv init -)"\nfi' >> ~/.bashrc
# restart shell
exec "$SHELL"
# install R v4.3.1
renv install 4.3.1
```

Clone this repository and install R packages (may take up to one hour):

```bash
# clone project repository
git clone https://github.com/boopsboops/ancistr.git
cd ancistr
# install R packages
Rscript -e "renv::restore()"
```


### USAGE

Run the analysis and generate the report (may take up to five minutes):

```bash
# creates file 'ancistr.html'
# previous analyses are kept in dated folders in 'results/'
quarto render ancistr.qmd
```


### CONTENTS (A-Z)

* **`assets/`** - Required files.
    - `ancistrus-coi.csv` - comma delimited table of in-house COI sequence data
    - `Ancistrus_sp(3).jpg` - photo of the bristlenose catfish
* **`renv/`** - Settings for the R environment.
* **`results/`** - Storage directory of previous analyses. Ignored by git.
* **`scripts/`** - R scripts.
    - `ancistr.R` - script to run the GenBank search and phylogenetic analysis
    - `load-libs-funs.R` - script loading up libraries and custom functions
* `ancistr.qmd` - quarto file to compile report and render as html
* `ancistr.html` - report output html file (ignored by git)
* `LICENSE` - Legal stuff
* `README.md` - This file
* `renv.lock` - R packages required and managed by renv
* `.gitignore` - files and directories ignored by git
* `.Rprofile` - activates renv
* `.R-version` - required version of R (used by renv-installer)
