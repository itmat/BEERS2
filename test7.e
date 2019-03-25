Traceback (most recent call last):
  File "bin/run_beers.py", line 34, in <module>
    args.func(args)
  File "/project/itmatlab/for_soum/BEERS2.0/beers/controller.py", line 78, in run_expression_pipeline
    self.input_samples)
  File "/project/itmatlab/for_soum/BEERS2.0/beers/expression/expression_pipeline.py", line 552, in main
    pipeline.execute()
  File "/project/itmatlab/for_soum/BEERS2.0/beers/expression/expression_pipeline.py", line 491, in execute
    outcome = beagle.execute(self.beagle_file_path, seeds["beagle"])
  File "/project/itmatlab/for_soum/BEERS2.0/beers/expression/beagle.py", line 20, in execute
    result = subprocess.call(command, shell=True)
  File "/opt/software/python/3.6.3/lib/python3.6/subprocess.py", line 269, in call
    return p.wait(timeout=timeout)
  File "/opt/software/python/3.6.3/lib/python3.6/subprocess.py", line 1457, in wait
    (pid, sts) = self._try_wait(0)
  File "/opt/software/python/3.6.3/lib/python3.6/subprocess.py", line 1404, in _try_wait
    (pid, sts) = os.waitpid(self.pid, wait_flags)
KeyboardInterrupt
