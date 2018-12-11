import os
import itertools
from collections import namedtuple
from beers.beers_exception import BeersException


Adapter = namedtuple('Adapter', ['label', 'end', 'sequence'])


class AdapterGenerator:

    next_adapter_labels_index = 0
    available_adapters = {}

    @staticmethod
    def generate_adapters(adapter_kit):
        if not AdapterGenerator.available_adapters:
            adapters = []
            current_file_path = os.path.abspath(os.path.dirname(__file__))
            adapter_file_path = os.path.join(current_file_path, f"../../resources/{adapter_kit}")
            if not os.path.exists(adapter_file_path):
                raise BeersException(f"No adapter kit resource is found at {adapter_file_path}.")
            with open(adapter_file_path) as adapter_file:
                for header, content in itertools.zip_longest(*[adapter_file] * 2):
                    fields = header[1:].split("-")
                    end = 3 if '3' in fields[1] else 5
                    label = fields[0][1:].rstrip()
                    sequence = content.rstrip()
                    adapter = Adapter(label, end, sequence)
                    adapters.append(adapter)
            labels_5 = [adapter.label for adapter in adapters if adapter.end == 5]
            labels_3 = [adapter.label for adapter in adapters if adapter.end == 3]
            labels = list(itertools.product(labels_5, labels_3))
            # Note - it is not expected that more than one adapter kit will be used in a given run.
            AdapterGenerator.available_adapters = {'adapters': adapters, 'labels': labels, 'current_index': 0}

    @staticmethod
    def get_unique_adapter_labels():
        adapter_data = AdapterGenerator.available_adapters
        adapter_labels = adapter_data['labels'][adapter_data['current_index']]
        adapter_data['current_index'] += 1
        if adapter_data['current_index'] >= len(adapter_data['labels']):
            raise BeersException("Too few adapter sequences available for the number of samples provided.")
        return adapter_labels

    @staticmethod
    def get_adapter_sequences_from_labels(adapter_labels):
        adapter_data = AdapterGenerator.available_adapters
        adapter_5_prime = [adapter.sequence for adapter in adapter_data['adapters']
                           if adapter.label == adapter_labels[0]][0]
        adapter_3_prime = [adapter.sequence for adapter in adapter_data['adapters']
                           if adapter.label == adapter_labels[1]][0]
        return adapter_5_prime, adapter_3_prime


if __name__ == '__main__':
    adapter_generator = AdapterGenerator("TruSeq_adapter_sequences_with_barcodes.MiSeq_HiSeq2000_HiSeq2500.fa")



