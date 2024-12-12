#####################################################################################################################
####               Google cloud script:  Convert genotype density file to bcf file on Google Cloud               ####
#####################################################################################################################

## 1. Convert .idat to .bcf using the MoChA pipeline on Google Cloud 
#to start the cromwell server on google VM
(java -XX:MaxRAMPercentage=90 -Dconfig.file=cromwell.conf -jar cromwell-85.jar server &)
java -jar cromwell-85.jar submit mocha.wdl -i pakistan_gsa_2433_vcf.json -o options_2433_gsa_vcf.json

# create the pakistan_gsa_2433_vcf.json
{
  "mocha.sample_set_id": "pakistan_gsa",
  "mocha.mode": "idat",
  "mocha.target": "vcf",
  "mocha.realign": true,
  "mocha.max_win_size_cm": 300.0,
  "mocha.overlap_size_cm": 5.0,
  "mocha.ref_name": "GRCh38",
  "mocha.ref_path": "gs://bigdeli-working/pakistan/GRCh38",
  "mocha.manifest_path": "gs://bigdeli-working/pakistan/manifests",
  "mocha.batch_tsv_file": "gs://bigdeli-working/pakistan/samples/pakistan_gsa.batch.tsv",
  "mocha.sample_tsv_file": "gs://bigdeli-working/pakistan/samples/pakistan_gsa.sample.tsv",
  "mocha.data_path": "gs://bigdeli-working/pakistan/idats"
}


# options_2433_gsa_vcf.json
{
  "delete_intermediate_output_files": true,
  "final_workflow_outputs_dir": "gs://bigdeli-working/pakistan/cromwell/outputs/gsa_2433_vcf",
  "use_relative_output_paths": true,
  "final_workflow_log_dir": "gs://bigdeli-working/pakistan/cromwell/wf_logs",
  "final_call_logs_dir": "gs://bigdeli-working/pakistan/cromwell/call_logs"
}
