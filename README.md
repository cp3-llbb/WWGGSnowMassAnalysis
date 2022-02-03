# WWGGSnowMassAnalysis

Repository for the HH->WWGG Snowmass analysis. 

`bambooRun -m CleanMerged.py:SnowmassExample YML/tautauGG_for_latex_table.yml -o OUTPUT_NAME` 

To test, one can use `--MaxFiles=1` to run with 1 root file.

To submit to slurm, only add `--distributed=driver`

`bambooRun --distributed=driver -m CleanMerged.py:SnowmassExample YML/tautauGG_for_latex_table.yml -o OUTPUT_NAME` 

If some jobs fail, resubmit with: 

`sbatch --array=JOB_NUMBER --export=ALL --licenses=cms_storage:3 path_to_the_output_folder/batch/slurmSubmission.sh`

and after those get completed:

`bambooRun --distributed=finalize -m CleanMerged.py:SnowmassExample YML/tautauGG_for_latex_table.yml -o OUTPUT_NAME `

## Options 

`--mvaSkim` to skim and make a TTree with skimmed variables
`--mvaEval` to import a DNN model in bamboo and evaluate on samples.
`--onlypost` if you already ran and have your outputs but willing to make changes to plots or YML details

