from abc import ABC, abstractmethod

class AbstractPipelineStep(ABC):
    """
    Abstract class defining the minimal methods required by a step in any pipeline
    of the BEERS suite of tools.
    """

    @abstractmethod
    def execute(self):
        """
        Entry point into the pipeline step.
        """
        pass

    @abstractmethod
    def validate(self):
        """
        Checks validity of parameters used to instantiate the pipeline step.

        Returns
        -------
        boolean
            True  - All parameters required to run this step were provided and
                    are within valid ranges.
            False - One or more of the paramters is missing or contains an invalid
                    value.
        """
        pass
