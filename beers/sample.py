class Sample:

    next_sample_id = 1

    def __init__(self, sample_id, sample_name, input_file_path, gender, adapter_labels):
        self.sample_id = sample_id
        self.sample_name = sample_name
        self.input_file_path = input_file_path
        self.gender = gender
        self.adapter_labels = adapter_labels

    def __str__(self):
        return f"sample id: {self.sample_id}, sample name: {self.sample_name},\n" \
               f"input_file_path: {self.input_file_path},\n" \
               f"gender: {self.gender}, adapter labels: {self.adapter_labels}"
