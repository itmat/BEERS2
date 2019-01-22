class Sample:
    """
    This object represents a physical sample in the BEER simulator run.  There can be many samples, the data for
    which may be spread across many molecule/cluster packets.  Each packet includes a pointer to the sample from which
    the packet data is derived.  Samples derive from the basic input of the overall pipeline and are created by the
    controller just prior to running the expression pipeline.  Since the sample is a component of both the molecule
    packet and the cluster packet, it must be serializable.
    """

    next_sample_id = 1

    def __init__(self, sample_id, sample_name, input_file_paths, adapter_sequences, gender=None):
        """
        The gender may not be known upon instantiation.
        :param sample_id: An integer uniquely indentifying the sample
        :param sample_name: The sample name, which is really the original sample file name
        :param input_file_path: The absolute input path to the sample.
        :param gender: 'male' or 'female'
        :param adapter_sequences: a list containing the 5' adapter sequence followed by the 3' adapter sequence.
        """
        # TODO Should adapter_sequences be replaced with the Adapter namedtuple, which does contain more information
        self.sample_id = sample_id
        self.sample_name = sample_name
        self.input_file_paths = input_file_paths
        self.gender = gender
        self.adapter_sequences = adapter_sequences

    def __str__(self):
        return f"sample id: {self.sample_id}, sample name: {self.sample_name},\n" \
               f"input_file_paths: {self.input_file_paths},\n" \
               f"gender: {self.gender}, adapter sequences: {self.adapter_sequences}"

    def serialize(self):
        """
        Single tab delimited string containing the sample attributes.  The adapter sequences are separated by a
        comma.
        :return: single string representing the sample data
        """
        adapter_sequences_tuple = ",".join([str(adapter_sequence) for adapter_sequence in self.adapter_sequences])
        return f"{self.sample_id}\t{self.sample_name}\t{self.input_file_paths}\t{self.gender}\t{adapter_sequences_tuple}"

    @staticmethod
    def deserialize(data):
        """
        Re-render the provided data back into a Sample object.  Any leading hash if stripped prior to unpacking the
        data string
        :param data: The serialized version of the sample object.
        :return: Sample object populated with the serialized data.
        """
        data = data[1:] if(data.startswith("#")) else data
        sample_id, sample_name, input_file_paths, gender, adapter_sequences_str = data.rstrip().split('\t')
        adapter_labels = tuple(str(adapter_label) for adapter_label in adapter_sequences_str.split(','))
        return Sample(sample_id, sample_name, input_file_paths, adapter_labels, gender)
