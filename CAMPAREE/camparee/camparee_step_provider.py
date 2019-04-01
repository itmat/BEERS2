from camparee.camparee_utils import CampareeException
from camparee.abstract_camparee_step import AbstractCampareeStep
#Add classes supporting additional pipeline steps here.
from camparee.variants_finder import VariantsFinderStep

class CampareeStepProvider:
    """
    Short summary.

    Attributes
    ----------
    __steps : dict
        Dictionary mapping pipeline step name, which is accessible and used by
        the rest of the code base, to the corresponding camparee_step. The
        camparee_step class must extend the AbstractCampareeStep class.

    """

    def __init__(self):
        self.__steps = {}

    def register_step(self, step_name, step_interface, package_name="camparee"):
        """
        Add interface to a given camparee step so it's accessible and useable by
        the rest of the code base.

        Parameters
        ----------
        step_name : string
            Name of the step corresponding to a specific pipeline step interface.
        step_interface : AbstractCampareeStep
            Class provding an interface to the camparee pipeline step.
        package_name : string
            Name of the package from which to load the interface class. [Default:
            camparee].

        """
        if issubclass(step_interface, AbstractCampareeStep):
            self.__steps[step_name] = step_interface
        else:
            raise CampareeStepProviderException(f'The {type(step_interface).__name__} class must'
                                                ' be a subclass of AbstractCampareeStep.')

    def list_supported_camparee_steps(self):
        """
        Return list of camparee_steps currently registered for use.

        Returns
        -------
        list
            Camparee_steps currently registered for use.

        """
        return list(self.__steps.keys())

    def get(self, step_name):
        """
        Return camparee_step class corresponding to the given step name.

        Parameters
        ----------
        step_name : string
            Name of the step corresponding to a specific camparee step interface.

        Returns
        -------
        AbstractCampareeStep
            Class providing interface to the camparee step.

        """
        step_interface = self.__steps.get(step_name)
        if not step_interface:
            raise CampareeStepProviderException(f'{step_name} is not a supported mode.\n'
                                                'Please select one of {",".join(self.__steps.keys())}.\n')
        return step_interface


class CampareeStepProviderException(CampareeException):
    pass


#As new step interfaces are implemented, they must be imported above and registered
#here, tying them to their associated step names.
CAMPAREE_STEPS = CampareeStepProvider()
CAMPAREE_STEPS.register_step("VariantsFinderStep", VariantsFinderStep)
