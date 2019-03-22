from beers.beers_exception import BeersException
from beers.abstract_job_scheduler import AbstractJobScheduler
#Add classes supporting additional job schedulers here.
from beers.lsf_job_scheduler import LsfJobScheduler

class JobSchedulerProvider:
    """
    Factory class responsible for providing the correct job scheduler interface
    given a dispatcher mode. This code generalization helps separate specific
    implementation details from the rest of the code base, and will make it
    easier to add support for other job schedulers in the future.

    Other code cann access the job_scheduler_provider framework using the
    following code:
        beers.job_scheduler_provider.SCHEDULERS.*

    Attributes
    ----------
    _schedulers : dict
        Dictionary mapping dispatcher_mode, which is accessible and used by the
        rest of the code base, to the corresponding job_scheduler class. The
        job_scheduler class must extend the AbstractJobScheduler class.

    """

    def __init__(self):
        self._schedulers = {}

    def register_scheduler(self, dispatcher_mode, scheduler):
        """
        Add interface to a given job scheduler so it's accessible and useable by
        the rest of the code base.

        Parameters
        ----------
        dispatcher_mode : string
            Mode corresponding to a specific job scheduler interface.
        scheduler : AbstractJobScheduler
            Class provding an interface to the job scheduler.

        """
        if issubclass(scheduler, AbstractJobScheduler):
            self._schedulers[dispatcher_mode] = scheduler
        else:
            raise JobSchedulerException(f'The {type(scheduler).__name__} class must'
                                        ' be a subclass of AbstractJobScheduler.')

    def list_supported_schedulers(self):
        """
        Return list of dispatcher_modes currently registered for use.

        Returns
        -------
        list
            Dispatcher_modes currently registered for use.

        """
        return list(self._schedulers.keys())

    #TODO: Change "dispatcher_mode" variable throught the software to more
    #accurately reflect what it represents. The dispatcher class was going to
    #handle submitting jobs, but has since been superseded by the more generic
    #job_scheduler framework. We should probably call it something like
    #"scheduler_mode" to better reflect where it's currently seeing the most use.
    #The dispatcher is still being used by the library prep pipeline and we'll
    #need to give it a thorough review before determining if/how we can remove it.
    def get(self, dispatcher_mode):
        """
        Return job_scheduler class corresponding to the given dispatcher mode.

        Parameters
        ----------
        dispatcher_mode : string
            Mode corresponding to a specific job scheduler interface.

        Returns
        -------
        AbstractJobScheduler
            Class providing interface to the job scheduler.

        """
        scheduler = self._schedulers.get(dispatcher_mode)
        if not scheduler:
            raise JobSchedulerException(f'{dispatcher_mode} is not a supported mode.\n'
                                        'Please select one of {",".join(self._schedulers.keys())}.\n')
        return scheduler

class JobSchedulerException(BeersException):
    pass

#As new job schedulers are implemented, they must be imported above and registered
#here, typing them to their associated dispatcher mode.
SCHEDULERS = JobSchedulerProvider()
SCHEDULERS.register_scheduler("lsf", LsfJobScheduler)
