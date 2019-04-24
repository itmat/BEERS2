from beers_utils.general_utils import BeersUtilsException
from beers_utils.abstract_job_scheduler import AbstractJobScheduler
#Add classes supporting additional job schedulers here.
from beers_utils.lsf_job_scheduler import LsfJobScheduler
from beers_utils.sge_job_scheduler import SgeJobScheduler

class JobSchedulerProvider:
    """
    Factory class responsible for providing the correct job scheduler interface
    given a scheduler mode. This code generalization helps separate specific
    implementation details from the rest of the code base, and will make it
    easier to add support for other job schedulers in the future.

    Other code cann access the job_scheduler_provider framework using the
    following code:
        beers_utils.job_scheduler_provider.SCHEDULERS.*

    Attributes
    ----------
    _schedulers : dict
        Dictionary mapping scheduler_mode, which is accessible and used by the
        rest of the code base, to the corresponding job_scheduler class. The
        job_scheduler class must extend the AbstractJobScheduler class.

    """

    def __init__(self):
        self._schedulers = {}

    def register_scheduler(self, scheduler_mode, scheduler):
        """
        Add interface to a given job scheduler so it's accessible and useable by
        the rest of the code base.

        Parameters
        ----------
        scheduler_mode : string
            Mode corresponding to a specific job scheduler interface.
        scheduler : AbstractJobScheduler
            Class provding an interface to the job scheduler.

        """
        if issubclass(scheduler, AbstractJobScheduler):
            self._schedulers[scheduler_mode] = scheduler
        else:
            raise JobSchedulerException(f'The {type(scheduler).__name__} class must'
                                        ' be a subclass of AbstractJobScheduler.')

    def list_supported_schedulers(self):
        """
        Return list of scheduler_modes currently registered for use.

        Returns
        -------
        list
            Scheduler_modes currently registered for use.

        """
        return list(self._schedulers.keys())

    def get(self, scheduler_mode):
        """
        Return job_scheduler class corresponding to the given scheduler mode.

        Parameters
        ----------
        scheduler_mode : string
            Mode corresponding to a specific job scheduler interface.

        Returns
        -------
        AbstractJobScheduler
            Class providing interface to the job scheduler.

        """
        scheduler = self._schedulers.get(scheduler_mode)
        if not scheduler:
            raise JobSchedulerException(f'{scheduler_mode} is not a supported mode.\n'
                                        'Please select one of {",".join(self._schedulers.keys())}.\n')
        return scheduler

class JobSchedulerException(BeersUtilsException):
    pass

#As new job schedulers are implemented, they must be imported above and registered
#here, typing them to their associated scheduler mode.
SCHEDULERS = JobSchedulerProvider()
SCHEDULERS.register_scheduler("lsf", LsfJobScheduler)
SCHEDULERS.register_scheduler("sge", SgeJobScheduler)
