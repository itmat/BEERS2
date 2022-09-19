from beers_utils.molecule import Molecule
import numpy as np
import pickle
from timeit import default_timer as timer
from beers.utilities.adapter_generator import AdapterGenerator
import beers_utils.cigar


class AdapterLigationStep:

    name = "Adapter Ligation Step"

    def __init__(self, step_log_file_path, parameters, global_config):
        self.log_filename = step_log_file_path
        self.parameters = parameters
        self.global_config = global_config
        print(f"{self.name} instantiated")

    def execute(self, molecule_packet, rng):
        print(f"{self.name} starting")
        sample = molecule_packet.sample
        # Adapters combine a fixed sequence (specified in 'resources' config)
        # with barcodes that are sample-specific
        i5_barcode = self.global_config['samples'][str(sample.sample_id)]['barcodes']['i5']
        i7_barcode = self.global_config['samples'][str(sample.sample_id)]['barcodes']['i7']
        adapter_5_prime = self.global_config['resources']['pre_i5_adapter'] + i5_barcode + self.global_config['resources']['post_i5_adapter']
        # NOTE: on the 3' end, an "A" gets ligated onto the sequence first. This is done as a discrete step
        # in the actual TruSeq protocol, but we do it as part of this step here.
        adapter_3_prime = "A" + self.global_config['resources']['pre_i7_adapter'] + i7_barcode + self.global_config['resources']['post_i7_adapter']

        # Ligate the adapters onto each molecule
        with open(self.log_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in molecule_packet.molecules:
                sequence = molecule.sequence
                molecule.sequence = adapter_5_prime + sequence + adapter_3_prime
                molecule.start = 1
                molecule.cigar = f"{len(adapter_5_prime)}S{len(sequence)}M{len(adapter_3_prime)}S"
                new_source_start, new_source_cigar, new_source_strand = beers_utils.cigar.chain(
                        molecule.start, molecule.cigar, "+",
                        molecule.source_start, molecule.source_cigar, molecule.source_strand
                )
                molecule.source_start = new_source_start
                molecule.source_cigar = new_source_cigar
                molecule.source_strand = new_source_strand
                log_file.write(molecule.log_entry())
        return molecule_packet

    @staticmethod
    def validate(parameters, global_config):
        errors = []
        for sample in global_config['samples'].values():
            if not 'i5' in sample['barcodes'] or not 'i7' in sample['barcodes']:
                errors.append(f"In sample {sample.sample_id}, expected to find i5 and i7 adpters in config.")

        adapters = ['pre_i5_adapter', 'pre_i7_adapter', 'post_i5_adapter', 'post_i7_adapter']
        if not all(adapter in global_config['resources'] for adapter in adapters):
            errors.append(f"In resources config, need to  specify all of {adapters}")

        return errors
