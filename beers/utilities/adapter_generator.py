import os
import itertools
from collections import namedtuple
from beers.beers_exception import BeersException


Adapter = namedtuple('Adapter', ['label', 'end', 'sequence'])

class AdapterGenerator:

    next_adapter_labels_index = 0

    def __init__(self, adapter_filename):
        self.adapters = []
        current_file_path = os.path.abspath(os.path.dirname(__file__))
        adapter_file_path = os.path.join(current_file_path, f"../../resources/{adapter_filename}")
        with open(adapter_file_path) as adapter_file:
            for header, content in itertools.zip_longest(*[adapter_file] * 2):
                    fields = header[1:].split("-")
                    end = 3 if '3' in fields[1] else 5
                    label = fields[0][1:].rstrip()
                    sequence = content.rstrip()
                    adapter = Adapter(label, end, sequence)
                    self.adapters.append(adapter)
        labels_5 = [adapter.label for adapter in self.adapters if adapter.end == 5]
        labels_3 = [adapter.label for adapter in self.adapters if adapter.end == 3]
        self.labels = list(itertools.product(labels_5, labels_3))

    def get_unique_adapter_labels(self):
        adapter_labels = self.labels[AdapterGenerator.next_adapter_labels_index]
        AdapterGenerator.next_adapter_labels_index += 1
        if AdapterGenerator.next_adapter_labels_index >= len(self.labels):
            raise BeersException("Too few adapter sequences available for the number of samples provided.")
        return adapter_labels


if __name__ == '__main__':
    adapter_generator = AdapterGenerator("TruSeq_adapter_sequences_with_barcodes.MiSeq_HiSeq2000_HiSeq2500.fa")



