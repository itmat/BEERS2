import shutil
import os

class GenomeAlignmentStep():

    name = "Genome Alignment Step"

    def __init__(self, logfile, data_directory_path, parameters):
        self.data_directory_path = data_directory_path

    def validate(self):
        return True

    def execute(self, sample, reference_genome):
        # TODO Stub!!  Hard-coded 2 alignment files to stub for STAR alignments.
        alignment_file_path = ''
        alignment_index_file_path = ''
        parent_dir = os.path.abspath(os.path.join(sample.input_file_path, os.pardir))
        if "1002" in sample.sample_name:
            alignment_file_path = os.path.join(parent_dir,"Test_data.1002_baseline.sorted.bam")
            alignment_index_file_path = os.path.join(parent_dir, "Test_data.1002_baseline.sorted.bai")
        if "1015" in sample.sample_name:
            alignment_file_path = os.path.join(parent_dir,"Test_data.1015_baseline.sorted.bam")
            alignment_index_file_path = os.path.join(parent_dir, "Test_data.1015_baseline.sorted.bai")
        shutil.copyfile(alignment_file_path,
                        os.path.join(self.data_directory_path, f'sample{sample.sample_id}', 'genome_alignment.bam'))
        shutil.copyfile(alignment_index_file_path,
                        os.path.join(self.data_directory_path, f'sample{sample.sample_id}', 'genome_alignment.bai'))



