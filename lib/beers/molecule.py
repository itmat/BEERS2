import re


class Molecule:

    next_molecule_id = 1 # Static variable for creating increasing molecule id's

    def __init__(self, molecule_id, sequence, start=None, cigar=None, strand="+"):
        self.molecule_id = molecule_id
        self.sequence = sequence
        self.start = start
        self.cigar = cigar
        self.strand = strand

    def poly_a_tail_length(self):
        match = re.search(r'(A+$)', self.sequence)
        return 0 if not match else len(match.group())

    def substitute(self, nucleotide, position):
        original_length = len(self.sequence)
        assert -1 <= position <= (original_length - 1), "Position must be along the molecule."
        self.sequence = self.sequence[:position] + nucleotide + self.sequence[position + 1:]
        self.start = None
        self.cigar = None

    def insert(self, insertion_sequence, position):
        # Position after which to insert
        # Use a position of -1 to prepend an insert
        original_length = len(self.sequence)
        insertion_length = len(insertion_sequence)
        assert -1 <= position <= original_length - 1, "Position must be along the molecule."
        if position == -1:
            self.cigar = f"{insertion_length}I{original_length}M"
            self.sequence = insertion_sequence + self.sequence
        elif position == original_length:
            self.cigar = f"{original_length}M{insertion_length}I"
            self.sequence = self.sequence + insertion_sequence
        else:
            lead_length = len(self.sequence[:position + 1])
            trail_length = original_length - position - 1
            self.cigar = f"{lead_length}M{insertion_length}I{trail_length}M"
            self.sequence = self.sequence[:position + 1] + insertion_sequence + self.sequence[position + 1:]

    def delete(self, deletion_length, position):
        # Position after which to delete
        # Use a position of -1 to delete from the 5' end
        original_length = len(self.sequence)
        assert -1 <= position <= original_length - 1, "Position must be along the molecule."
        assert original_length - deletion_length > 0, "Cannot delete the entire molecule in this way."
        if position == -1:
            self.cigar = f"{deletion_length}D{original_length}M"
            self.sequence = self.sequence[deletion_length:]
        elif position == original_length:
            self.cigar = f"{original_length}M{deletion_length}D"
            self.sequence = self.sequence[:original_length - deletion_length]
        else:
            lead_length = len(self.sequence[:position + 1])
            trail_length = original_length - deletion_length - lead_length
            self.cigar = f"{lead_length}M{deletion_length}D{trail_length}M"
            self.sequence = self.sequence[:position+1] + self.sequence[position + 1 + deletion_length:]

    def break_sequence(self, position, retain_3prime=True):
        # Position after which to break the molecule
        # For the present, assume that the 3 prime end is always the end retained.
        if self.strand == "+":
            self.start = self.start + position + 1
            self.sequence = self.sequence[position + 1:]
        else:
            self.sequence = self.sequence[:position + 1]
        self.cigar = f"{len(self.sequence)}M"

    def fragment(self, start,end):
        """ Return a smaller molecule from this molecule """
        assert start < end <= len(self.sequence)

        frag_sequence = self.sequence[start:end]
        frag_length = end - start
        frag_cigar = f"{frag_length}M" # Fragments match their parents
        frag_id = Molecule.new_id(self.molecule_id)

        frag = Molecule(frag_id, frag_sequence, start=start, cigar=frag_cigar, strand=self.strand)

        return frag

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return(
            str({"id": self.molecule_id,
                 "sequence": self.sequence,
                 "start": self.start,
                 "cigar": self.cigar,
                 "strand": self.strand})
        )

    def log_entry(self):
        return str(self.molecule_id) + "," + \
               str(self.sequence) + "," + str(self.start or '') + "," + str(self.cigar or '') + "," + self.strand

    @staticmethod
    def new_id(parent_id=""):
        new_id = Molecule.next_molecule_id
        Molecule.next_molecule_id += 1
        return f"{parent_id}.{new_id}"


if __name__ == "__main__":
    source = "AGTTCAAGCTTGCACTCTAG"
    source_length = len(source)
    molecule = Molecule("1", source, start="1", cigar=f"{source_length}M")
    print(molecule)
    molecule.substitute("C", 15)
    print(molecule)
    molecule.insert("AGT", 4)
    print(molecule)
    molecule.delete(2, 2)
    print(molecule)
