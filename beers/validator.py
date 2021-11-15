from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline
from beers.sequence.sequence_pipeline import SequencePipeline

class Validator:

    @staticmethod
    def validate(stage_names, configuration):
        for stage_name in stage_names:
            if stage_name == 'library_prep_pipeline':
                LibraryPrepPipeline.validate(configuration[stage_name], configuration)
            if stage_name == 'sequence_pipeline':
                SequencePipeline.validate(configuration[stage_name], configuration)


