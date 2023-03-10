import json
import numpy as np
import os

from typing import Annotated, Literal, Union
from pydantic import BaseModel, Field, conint, confloat, constr, validator

########################
#  General parameters  #
########################

Oligonucleotide = constr(regex=r"^(A|C|G|T)+$", strip_whitespace=True, to_upper=True)


class OutputConfiguration(BaseModel):
    output_fastq: bool = True
    output_sam: bool = True
    output_bam: bool = False
    full_logs: bool = False


class Barcodes(BaseModel):
    i5: Oligonucleotide
    i7: Oligonucleotide


class Sample(BaseModel):
    barcodes: Barcodes


class Resources(BaseModel):
    pre_i5_adapter: Oligonucleotide
    post_i5_adapter: Oligonucleotide
    pre_i7_adapter: Oligonucleotide
    post_i7_adapter: Oligonucleotide
    reference_genome_fasta: Path


class MoleculeMakerParameters(BaseModel):
    min_polyA_tail_length: conint(ge=0) = 40
    max_polyA_tail_length: conint(ge=0)


class GlobalConfiguration(BaseModel):
    samples: dict[constr(regex=r"^[0-9]+$"), Sample]
    molecule_maker_parameters: MoleculeMakerParameters
    resources: Resources

    @validator("molecule_maker_parameters")
    def check_bounds_valid(cls, values):
        assert (
            values.min_polyA_tail_length <= values.max_polyA_tail_length
        ), "ensure min_polyA_tail_length <= max_polyA_tail_length"
        return values

    @validator("samples")
    def check_barcode_lengths_agree(cls, samples):
        i5_barcode_lengths = [len(sample.barcodes.i5) for sample in samples.values()]
        i7_barcode_lengths = [len(sample.barcodes.i7) for sample in samples.values()]

        assert (
            len(set(i5_barcode_lengths)) == 1
        ), "ensure i5 barcodes are of the same length"

        assert (
            len(set(i7_barcode_lengths)) == 1
        ), "ensure i7 barcodes are of the same length"

        return samples


#########################
#  Library preparation  #
#########################


class MoleculeGenerationFromDistributionParameters(BaseModel):
    num_packets: conint(ge=0)
    num_molecules_per_packet: conint(ge=1)
    sample_data_directory: Path


class LibraryPreparationInput(BaseModel):
    directory_path: Path | None = None
    from_distribution_data: dict[
        constr(regex=r"^[0-9]+$"), MoleculeGenerationFromDistributionParameters
    ]


class PolyAStepParameters(BaseModel):
    breakpoint_prob_per_base: confloat(ge=0, le=1) = 0
    min_retention_prob: confloat(ge=0, le=1) = 0
    max_retention_prob: confloat(ge=0, le=1) = 1
    length_retention_prob: confloat(ge=0, le=1) = 1
    min_polya_tail_length: conint(ge=0) = 40


class PolyAStepConfiguration(BaseModel):
    step_name: Literal["polya_step.PolyAStep"]
    parameters: PolyAStepParameters

    @validator("parameters")
    def check_bound_constraints(cls, parameters):
        min_retention_prob = parameters.min_retention_prob
        max_retention_prob = parameters.max_retention_prob

        assert (
            min_retention_prob <= max_retention_prob
        ), "ensure min_retention_prob <= max_retention_prob"
        return parameters


class RiboZeroStepParameters(BaseModel):
    max_degrade_chance: confloat(ge=0, le=1)
    degrade_exponent: confloat(gt=0)
    degrade_entire_molecule: bool


class RiboZeroStepConfiguration(BaseModel):
    step_name: Literal["ribozero_step.RiboZeroStep"]
    parameters: RiboZeroStepParameters


class UniformFragmentMethod(BaseModel):
    method: Literal["uniform"]
    rate: confloat(gt=0)
    runtime: confloat(gt=0)
    min_frag_size: conint(gt=0)


class BetaFragmentMethod(BaseModel):
    method: Literal["beta"]
    rate: confloat(gt=0)
    runtime: confloat(gt=0)
    min_frag_size: conint(gt=0)
    beta_A: confloat(gt=0)
    beta_B: confloat(gt=0)
    beta_N: confloat(ge=0)


class FragmentStepConfiguration(BaseModel):
    step_name: Literal["fragment_step.FragmentStep"]
    parameters: UniformFragmentMethod | BetaFragmentMethod = Field(
        ..., discriminator="method"
    )


class PositionProbabilityMatrix(BaseModel):
    A: list[confloat(gt=0, lt=1)]
    C: list[confloat(gt=0, lt=1)]
    G: list[confloat(gt=0, lt=1)]
    T: list[confloat(gt=0, lt=1)]


class StrandSynthesisStepParameters(BaseModel):
    perfect_priming: bool
    position_probability_matrix: PositionProbabilityMatrix
    primes_per_kb: confloat(gt=0)

    @validator("position_probability_matrix")
    def check_probability_matrix_constraints(cls, p):
        assert (
            len(p.A) == len(p.C) == len(p.G) == len(p.T)
        ), "ensure probability lists for A, C, G, and T have equal length"

        assert all(
            np.isclose(sum(p_position), 1.0) for p_position in zip(*dict(p).values())
        ), f"ensure probabilities for A, C, G, and T at every position sum to 1"

        return p


class FirstStrandSynthesisStepConfiguration(BaseModel):
    step_name: Literal["first_strand_synthesis_step.FirstStrandSynthesisStep"]
    parameters: StrandSynthesisStepParameters


class SecondStrandSynthesisStepConfiguration(BaseModel):
    step_name: Literal["second_strand_synthesis_step.SecondStrandSynthesisStep"]
    parameters: StrandSynthesisStepParameters


class SizingStepParameters(BaseModel):
    min_length: conint(ge=0)
    max_length: conint(ge=0)
    select_all_start_length: conint(ge=0)
    select_all_end_length: conint(ge=0)

    @validator("select_all_start_length", always=True)
    def set_select_all_start_length_value(cls, select_all_start_length, values):
        return select_all_start_length or values["min_length"]

    @validator("select_all_end_length", always=True)
    def set_select_all_end_length_value(cls, select_all_end_length, values):
        return select_all_end_length or values["max_length"]


class SizingStepConfiguration(BaseModel):
    step_name: Literal["sizing_step.SizingStep"]
    parameters: SizingStepParameters

    @validator("parameters")
    def check_parameters(cls, parameters):
        assert (
            parameters.min_length
            <= parameters.select_all_start_length
            <= parameters.select_all_end_length
            <= parameters.max_length
        ), "ensure min_length <= select_all_start_length <= select_all_end_length <= max_length"

        return parameters


class AdapterLigationStepParameters(BaseModel):
    pass


class AdapterLigationStepConfiguration(BaseModel):
    step_name: Literal["adapter_ligation_step.AdapterLigationStep"]
    parameters: AdapterLigationStepParameters


class PCRAmplificationStepParameters(BaseModel):
    number_cycles: conint(ge=0, le=16)
    retention_percentage: confloat(gt=0, le=100)
    gc_bias_constant: float
    gc_bias_linear: float
    gc_bias_quadratic: float
    deletion_rate: confloat(ge=0, le=1)
    insertion_rate: confloat(ge=0, le=1)
    substitution_rate: confloat(ge=0, le=1)


class PCRAmplificationStepConfiguration(BaseModel):
    step_name: Literal["pcr_amplification_step.PCRAmplificationStep"]
    parameters: PCRAmplificationStepParameters


LibraryPreparationStepConfiguration = Annotated[
    Union[
        PolyAStepConfiguration,
        RiboZeroStepConfiguration,
        FragmentStepConfiguration,
        FirstStrandSynthesisStepConfiguration,
        SecondStrandSynthesisStepConfiguration,
        SizingStepConfiguration,
        AdapterLigationStepConfiguration,
        PCRAmplificationStepConfiguration,
    ],
    Field(discriminator="step_name"),
]


class LibraryPreparationPipelineConfiguration(BaseModel):
    input: LibraryPreparationInput
    steps: list[LibraryPreparationStepConfiguration]


################
#  Sequencing  #
################


class BridgeAmplificationStepParameters(BaseModel):
    cycles: conint(gt=0)
    substitution_rate: confloat(ge=0, le=1) = 0


class BridgeAmplificationStepConfiguration(BaseModel):
    step_name: Literal["bridge_amplification_step.BridgeAmplificationStep"]
    parameters: BridgeAmplificationStepParameters


class SequenceBySynthesisStepParameters(BaseModel):
    forward_is_5_prime: bool = True
    paired_ends: bool = True
    read_length: conint(gt=0)
    skip_rate: confloat(ge=0)
    drop_rate: confloat(ge=0)


class SequenceBySynthesisStepConfiguration(BaseModel):
    step_name: Literal["sequence_by_synthesis_step.SequenceBySynthesisStep"]
    parameters: SequenceBySynthesisStepParameters


class StepConfiguration(BaseModel):
    step_name: str
    parameters: dict


SequenceStepConfiguration = Annotated[
    Union[
        BridgeAmplificationStepConfiguration,
        SequenceBySynthesisStepConfiguration,
    ],
    Field(discriminator="step_name"),
]


class FlowcellGeometry(BaseModel):
    min_lane: int = 0
    max_lane: int = 10000
    min_tile: int = 0
    max_tile: int = 10000
    min_x: int = 0
    max_x: int = 10000
    min_y: int = 0
    max_y: int = 10000


class Flowcell(BaseModel):
    coordinate_strategy: Literal["random", "random_distinct"]
    flowcell_geometry: FlowcellGeometry
    lanes_to_use: list[int] | None

    @validator("flowcell_geometry")
    def check_flowcell_geometry(cls, flowcell):
        assert flowcell.min_lane <= flowcell.max_lane, "ensure min_lane <= max_lane"
        assert flowcell.min_tile <= flowcell.max_tile, "ensure min_tile <= max_tile"
        assert flowcell.min_x <= flowcell.max_x, "ensure min_x <= max_x"
        assert flowcell.min_y <= flowcell.max_y, "ensure min_y <= max_y"

        return flowcell

    @validator("lanes_to_use", always=True)
    def check_selected_lanes(cls, lanes, values):
        min_lane = values["flowcell_geometry"].min_lane
        max_lane = values["flowcell_geometry"].max_lane

        if not lanes:
            return list(range(min_lane, max_lane + 1))
        else:
            assert all(
                min_lane <= lane <= max_lane for lane in lanes
            ), "ensure all selected lanes are within [min_lane, max_lane] range"

            return lanes


class SequencePipelineConfiguration(BaseModel):
    steps: list[SequenceStepConfiguration]
    flowcell: Flowcell


########################
#  Full configuration  #
########################


class Configuration(BaseModel):
    seed: int
    output: OutputConfiguration
    library_prep_pipeline: LibraryPreparationPipelineConfiguration
    sequence_pipeline: SequencePipelineConfiguration
    global_config: GlobalConfiguration

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.check_all_samples_are_specified()

    def check_all_samples_are_specified(self):
        samples = self.global_config.samples.keys()
        missing_samples = set(
            self.library_prep_pipeline.input.from_distribution_data.keys()
        ) - set(samples)

        if missing_samples:
            location = "library_prep_pipeline -> input -> from_distribution_data"
            message = f"  ensure all samples listed are also specified in global_config - sample {missing_samples.pop()} missing in global_config"
            raise AssertionError(location + "\n" + message)
