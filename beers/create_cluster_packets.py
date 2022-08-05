        stage_name = "sequence_pipeline"
        if setup:
            self.perform_setup(args, [self.controller_name, stage_name])
        input_directory_path = self.configuration[stage_name]["input"]["directory_path"]
        intermediate_directory_path = os.path.join(self.output_directory_path, self.controller_name)
        data_directory_path = os.path.join(intermediate_directory_path, CONSTANTS.DATA_DIRECTORY_NAME)
        cluster_packet_file_paths = []
        molecule_packet_file_paths = glob.glob(f'{input_directory_path}{os.sep}**{os.sep}*.txt', recursive=True)
        if not molecule_packet_file_paths:
            raise ControllerValidationException(f"No molecule packet files were found in the given input directory, "
                                                f"{input_directory_path}, or any of its subdirectories")
        file_count = len(molecule_packet_file_paths)
        directory_structure = GeneralUtils.create_subdirectories(file_count, data_directory_path)
        self.setup_flowcell()
        for molecule_packet_filename in molecule_packet_file_paths:
            print(f"Loading {molecule_packet_filename} to flowcell")
            molecule_packet = MoleculePacket.deserialize(molecule_packet_filename)
            cluster_packet = self.flowcell.load_flowcell(molecule_packet)
            cluster_packet_filename = f"cluster_packet_start_pkt{cluster_packet.cluster_packet_id}.gzip"
            subdirectory_list = \
                GeneralUtils.get_output_subdirectories(cluster_packet.cluster_packet_id, directory_structure)
            data_subdirectory_path = os.path.join(data_directory_path, *subdirectory_list)
            cluster_packet_file_path = os.path.join(data_subdirectory_path, cluster_packet_filename)
            cluster_packet.serialize(cluster_packet_file_path)
            cluster_packet_file_paths.append(cluster_packet_file_path)
            print(f"Intermediate loaded - process RAM at {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1E6} GB")

