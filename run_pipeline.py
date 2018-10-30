from beers.library_prep.pipeline import Pipeline

lib_prep = Pipeline(config_filename="config/config.json")

lib_prep.validate()
lib_prep.execute()
