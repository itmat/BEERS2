import shutil
import os

class GenomeAlignmentStep():

    name = "Genome Alignment Step"

    sample_map = {
        "1002": "Test_data.1002_baseline.sorted.bam",
        "1015": "Test_data.1015_baseline.sorted.bam",
        "ILB_9577": "ILB_9577.Aligned.sortedByCoord.out.bam",
        "ILB_9578": "ILB_9578.Aligned.sortedByCoord.out.bam",
        "ILB_9583": "ILB_9583.Aligned.sortedByCoord.out.bam",
        "UNT_9576": "UNT_9576.Aligned.sortedByCoord.out.bam",
        "UNT_9580": "UNT_9580.Aligned.sortedByCoord.out.bam",
        "UNT_9584": "UNT_9584.Aligned.sortedByCoord.out.bam",
        "9576": "aligned/baby_sample1.bam",
        "9577": "aligned/baby_sample2.bam",
        "9578": "aligned/baby_sample3.bam"
    }


    def __init__(self, logfile, data_directory_path, parameters):
        self.data_directory_path = data_directory_path

    def validate(self):
        return True

    def execute(self, sample, reference_genome):
        # TODO Stub!!  Hard-coded alignment files to stub for STAR alignments.
        alignment_file_path = ''
        alignment_index_file_path = ''
        parent_dir = os.path.abspath(os.path.join(sample.input_file_paths[0], os.pardir))
        for key, value in GenomeAlignmentStep.sample_map.items():
            print(sample.sample_name)
            if key in sample.sample_name:
                alignment_file_path = os.path.join(parent_dir, value)
                alignment_index_file_path = os.path.join(parent_dir, value + ".bai")
                break
        shutil.copyfile(alignment_file_path,
                        os.path.join(self.data_directory_path, f'sample{sample.sample_id}', 'genome_alignment.bam'))
        shutil.copyfile(alignment_index_file_path,
                        os.path.join(self.data_directory_path, f'sample{sample.sample_id}', 'genome_alignment.bai'))



