from beers_utils.molecule import Molecule
from beers_utils.sample import Sample
from beers_utils.molecule_packet import MoleculePacket
from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline
from beers.library_prep.polya_step import PolyAStep
from beers.library_prep.first_strand_synthesis_step import FirstStrandSynthesisStep
from beers.library_prep.second_strand_synthesis_step import SecondStrandSynthesisStep
from beers.library_prep.fragment_step import FragmentStep
from beers.library_prep.sizing_step import SizingStep
from beers.library_prep.adapter_ligation_step import AdapterLigationStep
from beers.library_prep.pcr_amplification_step import PCRAmplificationStep
import numpy

def make_molecule_packet(count, length, rng):
    def random_sequence(length):
        return ''.join(rng.choice(['A', 'C', 'G', 'T'], size=length))

    return MoleculePacket(
        molecules = [
            Molecule(
                molecule_id = "5",
                sequence = random_sequence(length),
                start = 5,
                cigar = f"{length}M",
                strand = '+',
                transcript_id = "ENSMUSG0000000000",
                source_start = rng.integers(1, 5_000_000_000),
                source_cigar =  f"{length}M",
                source_strand = rng.choice(["-", "+"]),
                source_chrom = "X",
            )
            for i in range(10)
        ],
        molecule_packet_id = 5,
        sample = Sample("1", "sample1", None, None, None)
    )



def test_PolyA(tmp_path):
    # Setup
    rng = numpy.random.default_rng(0)

    molecule_packet = make_molecule_packet(count = 10, length= 3_000, rng = rng)
    for molecule in molecule_packet.molecules[:5]:
        # With polyA tails for the first 5
        molecule.sequence = molecule.sequence + "A"*500
        molecule.source_cigar = molecule.source_cigar + "500S"
    original_molecules = molecule_packet.molecules

    step = PolyAStep(
        step_log_file_path = tmp_path / "log.txt",
        parameters = {
            'breakpoint_prob_per_base': 0.01,
            'max_retention_prob': 1.0,
            'min_retention_prob': 0.0,
            'min_polya_tail_length': 40,
            'length_retention_prob': 0.05,
        },
        global_config = {}
    )

    # Run the step
    output = step.execute(molecule_packet, rng)

    # Exactly 5 had PolyA tails
    assert len(output.molecules) == 5
    for molecule, original in zip(output.molecules, original_molecules):
        # And all should be a subset of the earlier molecule
        assert len(molecule.sequence) <= len(original.sequence)
        assert molecule.sequence in original.sequence

def test_FirstStrandSynthesis(tmp_path):
    # Setup
    rng = numpy.random.default_rng(0)

    molecule_packet = make_molecule_packet(count = 10, length= 3_000, rng = rng)
    original_molecules = molecule_packet.molecules

    step = FirstStrandSynthesisStep(
        log_file = tmp_path / "log.txt",
        parameters = {
            "perfect_priming": False,
            "position_probability_matrix": {
                'A': [0.50, 0.1, 0.40, 0.30, 0.25, 0.15],
                'C': [0.20, 0.5, 0.3 , 0.25, 0.25, 0.15],
                'G': [0.15, 0.1, 0.15, 0.25, 0.25, 0.20],
                'T': [0.15, 0.3, 0.15, 0.20, 0.25, 0.50],
            },
            "primes_per_kb": 50
        },
        global_config = {}
    )

    # Run the step
    output = step.execute(molecule_packet, rng)

    # Retain all molecules
    assert len(output.molecules) == 10

    for molecule, original in zip(output.molecules, original_molecules):
        # And all should be a subset of the original molecule
        # but reverse complemented
        assert len(molecule.sequence) <= len(original.sequence)


def test_SecondStrandSynthesis(tmp_path):
    # Setup
    rng = numpy.random.default_rng(0)

    molecule_packet = make_molecule_packet(count = 10, length= 3_000, rng = rng)
    original_molecules = molecule_packet.molecules

    step = SecondStrandSynthesisStep(
        log_file = tmp_path / "log.txt",
        parameters = {
            "perfect_priming": False,
            "position_probability_matrix": {
                'A': [0.50, 0.1, 0.40, 0.30, 0.25, 0.15],
                'C': [0.20, 0.5, 0.3 , 0.25, 0.25, 0.15],
                'G': [0.15, 0.1, 0.15, 0.25, 0.25, 0.20],
                'T': [0.15, 0.3, 0.15, 0.20, 0.25, 0.50],
            },
            "primes_per_kb": 50
        },
        global_config = {}
    )

    # Run the step
    output = step.execute(molecule_packet, rng)

    # Retain all molecules
    assert len(output.molecules) == 10

    for molecule, original in zip(output.molecules, original_molecules):
        # And all should be a subset of the original molecule
        # but reverse complemented
        assert len(molecule.sequence) <= len(original.sequence)


def test_FragmentStep(tmp_path):
    # Setup
    rng = numpy.random.default_rng(0)

    molecule_packet = make_molecule_packet(count = 10, length= 3_000, rng = rng)
    original_molecules = molecule_packet.molecules

    step = FragmentStep(
        logfile = tmp_path / "log.txt",
        parameters = {
            "method": "uniform",
            "lambda": 0.005,
            "runtime": 1,
            "min_frag_size": 20,
        },
        global_config = {}
    )

    # Run the step
    output = step.execute(molecule_packet, rng)

    for molecule in output.molecules:
        # Original molecules were length 3_000
        # now we should have a subset
        assert len(molecule) <= 3_000
        # Each fragment is a subste of one of the starting molecules
        assert any(molecule.sequence in original_molecule.sequence
                        for original_molecule in original_molecules)
    # And the total length should be no bigger than the starting length
    total_frag_bases = sum(len(mol.sequence) for mol in output.molecules)
    total_original_bases = sum(len(mol.sequence) for mol in original_molecules)
    assert total_frag_bases <= total_original_bases

def test_SizingStep(tmp_path):
    # Setup
    rng = numpy.random.default_rng(0)

    small = make_molecule_packet(count = 10, length= 50, rng = rng)
    medium = make_molecule_packet(count = 10, length= 250, rng = rng)
    large = make_molecule_packet(count = 10, length= 1_000, rng = rng)
    molecule_packet = MoleculePacket(
        molecules = small.molecules + medium.molecules + large.molecules,
        sample = None,
        molecule_packet_id = 5
    )

    step = SizingStep(
        step_log_file_path = tmp_path / "log.txt",
        parameters = {
            "min_length": 100,
            "max_length": 400,
            "select_all_start_length": 200,
            "select_all_end_length": 300,
        },
        global_config = {}
    )

    # Run the step
    output = step.execute(molecule_packet, rng)

    # Retain exactly the medium length ones
    assert set(output.molecules) == set(medium.molecules)

def test_AdapterLigationStep(tmp_path):
    # Setup
    rng = numpy.random.default_rng(0)

    molecule_packet = make_molecule_packet(count = 10, length= 300, rng = rng)
    original_molecules = molecule_packet.molecules

    step = AdapterLigationStep(
        step_log_file_path = tmp_path / "log.txt",
        parameters = {
        },
        global_config = {
            "samples": {
                '1': dict( barcodes = dict( i5 = "AGCGCTAG", i7 = "AACCGCGG") ),
                '2': dict( barcodes = dict( i5 = "GATATCGA", i7 = "TTATAACC") ),
            },
            "resources": {
                "pre_i5_adapter": "AATGATACGGCGACCACCGAGATCTACAC",
                "post_i5_adapter": "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
                "pre_i7_adapter": "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
                "post_i7_adapter": "ATCTCGTATGCCGTCTTCTGCTTG",
            }
        }
    )

    # Run the step
    output = step.execute(molecule_packet, rng)

    assert len(output.molecules) == len(original_molecules)
    for molecule, original in zip(output.molecules, original_molecules):
        # And all should be a subset of the earlier molecule
        assert len(molecule.sequence) <= len(original.sequence)
        assert molecule.sequence in original.sequence



def test_PCRAmplificationStep(tmp_path):
    # Setup
    rng = numpy.random.default_rng(0)

    molecule_packet = make_molecule_packet(count = 10, length= 300, rng = rng)
    original_molecules = molecule_packet.molecules

    step = PCRAmplificationStep(
        step_log_file_path = tmp_path / "log.txt",
        parameters = {
            "number_cycles": 10,
            "retention_percentage": 0.08,
            "gc_bias_constant": 1.0,
            "gc_bias_linear": 0.0,
            "gc_bias_quadratic": 0.0,
            "deletion_rate": 0.0001,
            "insertion_rate": 0.0001,
            "substitution_rate": 0.001,
        },
        global_config = {
        }
    )

    # Run the step
    output = step.execute(molecule_packet, rng)

    # TODO: not sure what to verify other than that it ran
    #       PCR can change sequence, lengths, and either increase or decrease
    #       the number of molecules, all depending upon random chance.
