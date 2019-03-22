Traceback (most recent call last):
  File "bin/run_beers.py", line 34, in <module>
    args.func(args)
  File "/project/itmatlab/for_soum/BEERS2.0/beers/controller.py", line 78, in run_expression_pipeline
    self.input_samples)
  File "/project/itmatlab/for_soum/BEERS2.0/beers/expression/expression_pipeline.py", line 552, in main
    pipeline.execute()
  File "/project/itmatlab/for_soum/BEERS2.0/beers/expression/expression_pipeline.py", line 229, in execute
    self.star_file_path, self.dispatcher_mode)
  File "/project/itmatlab/for_soum/BEERS2.0/beers/expression/genome_alignment.py", line 51, in execute
    read_files = ' '.join(sample.input_file_paths)
AttributeError: 'Sample' object has no attribute 'input_file_paths'
