# WWGGSnowMassAnalysis

Repository for the HH->WWGG Snowmass analysis. 

The analysis is being performed with [bamboo framework](https://bamboo-hep.readthedocs.io/en/latest/).

`bambooRun -m CleanMerged.py:SnowmassExample YML/tautauGG_for_latex_table.yml -o OUTPUT_NAME` 

To test, one can use `--MaxFiles=1` to run with 1 root file.

To submit to slurm, only add `--distributed=driver`

`bambooRun --distributed=driver -m CleanMerged.py:SnowmassExample YML/tautauGG_for_latex_table.yml -o OUTPUT_NAME` 

If some jobs fail, resubmit with: 

`sbatch --array=JOB_NUMBER --export=ALL --licenses=cms_storage:3 path_to_the_output_folder/batch/slurmSubmission.sh`

and after those get completed:

`bambooRun --distributed=finalize -m CleanMerged.py:SnowmassExample YML/tautauGG_for_latex_table.yml -o OUTPUT_NAME `

A tiny problem that we experienced is that when you want to read files from eos spaces or some other storage, there was a problem with the proxy transfer. To avoid the `Could not open file` error, do:

```
cp $X509_USER_PROXY ~/.x509_proxy
export X509_USER_PROXY=~/.x509_proxy
```
So that the generated proxy can be transferred to the machines where your jobs are submitted to. 

## Options 

`--mvaSkim` to skim and make a TTree with skimmed variables.

`--mvaEval` to import a DNN model in bamboo and evaluate on samples.

`--onlypost` if you already ran and have your outputs but willing to make changes to plots or YML details

