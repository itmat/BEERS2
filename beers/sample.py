class Sample:

    next_sample_id = 1

    def __init__(self, sample_id, sample_name, input_file_path, gender):
        self.sample_id = sample_id
        self.sample_name = sample_name
        self.input_file_path = input_file_path
        self.gender = gender
        self.adapter_5_prime = ''
        self.adapter_3_prime = ''

    def add_adapters(self, adapter_5_prime, adapter_3_prime):
        self.adapter_5_prime = adapter_5_prime
        self.adapter_3_prime = adapter_3_prime
