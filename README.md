# thumper


Thumper is a snakemake workflow that uses [sourmash](https://github.com/sourmash-bio/sourmash) for genome classification.
It primarily leverages sourmash gather (draft manuscript here(https://github.com/dib-lab/2020-paper-sourmash-gather)) and
the `sourmash tax` subcommand introduced in `sourmash 4.2.1`.

For microbial genomes, we have been exploring the utility of protein k-mers for classification.

## Development Install

We suggest installing in an isolated conda environment.

To test and develop, first install Miniconda and clone this repo to your system. Then install an editable version of thumper.

```
git clone https://github.com/bluegenes/thumper.git
cd thumper
conda env create -f environment.yml
conda activate thumper
pip install -e "."
```
> Note: [`mamba`](https://github.com/mamba-org/mamba) works well here too!

## Run a demo

```
thumper run demo/protein.demo.yaml genome_classify
```

### Name Origin
Thumper is part of the `sourmash` family of whiskey-adjacent software names.

Thumper: One of the types of stills used to accomplish the second distillation of American whiskey. It effectively removes impurities and concentrates the alcohol even further. “Low wines” go in; “high wines” come out. Thumpers differ from doublers in that the low wines enter a thumper as vapors that are bubbled through water, causing the stills to make a thumping sound; a doubler makes no distinctive noise since the low wines enter in condensed, liquid form.
