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

    def serialize(self):
       adapter_labels_tuple = ",".join([str(adapter_label) for adapter_label in self.adapter_labels])
       return f"{self.sample_id}\t{self.sample_name}\t{self.input_file_path}\t{self.gender}\t{adapter_labels_tuple}"

    @staticmethod
    def deserialize(data):
        data = data[1:] if(data.startswith("#")) else data
        sample_id, sample_name, input_file_path, gender, adapter_labels_str = data.rstrip().split('\t')
        adapter_labels = tuple(str(adapter_label) for adapter_label in adapter_labels_str.split(','))
        return Sample(sample_id, sample_name, input_file_path, gender, adapter_labels)
