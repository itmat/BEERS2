#!/usr/bin/python

import json
import importlib

class Pipeline:

    def __init__(self):
        with open("../../config/config.json", "r+") as config_file:
            config_json = json.load(config_file)
            self.steps = []
            for step_name in config_json["pipeline_steps"]:
                module = importlib.import_module(step_name)
                step_class = getattr(module, step_name)
                self.steps.append(step_class())


    def validate(self, **kwargs):
        if not all([step.validate() for step in self.steps]):
            raise BeersValidationException("Validation error: see log for details.")

    def execute(self):
        sample = "Original RNA input with many molecule objects in a list"
        for step in self.steps:
            sample = step.execute(sample)

    @staticmethod
    def main():
        pipeline = Pipeline()

        pipeline.validate()
        pipeline.execute()


class BeersValidationException(Exception):
    pass


if __name__ == "__main__":
    Pipeline.main()
