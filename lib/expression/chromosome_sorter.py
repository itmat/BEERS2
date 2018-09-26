import roman
import sys
import re
import argparse


class ChromosomeSort:
    """
    Provides a facility for sorting a list of chromosome names by means of a ChromosomeName class.
    The ordering is as follows:
    0.  no character string (a test needed because of recursion - see below)
    1.  leading digits (sorted numerically)
    2.  leading roman numerals (sorted by arabic number equivalents)
    3.  leading names potentially identified as gender related chromosomes (sorted x,y,m)
    4.  leading 'chr'
    5.  leading alphabetic characters (sorted by dictionary order)
    6.  leading non-alphanumeric characters (sorted by ascii order)

    Case is disreguarded.  For any of the cases above, in the case of equivalence of these leading characters, any
    following character strings that start with a character different from the leading character classification,
    are used to create new ChromosomeName classes and then compared - essentially recusing through the string
    representing the ChromosomeName.
    """

    def __init__(self, chromosome_name_filename):
        """
        Holds the file containing the list of chromosome names.
        :param chromosome_name_filename:
        """
        self.chromosome_name_filename = chromosome_name_filename

    def sort_chromosome_names(self):
        """
        Sorts the chromosome names provided according to the rules stated in the class documentation.  Both the
        unsorted and sorted list are printed.
        """

        # Populate a list with the chromosome names (1 per line) from the provided file.  Each chromosome name is
        # used to instantiate a ChromosomeName object, which is capable or making comparisons.
        chromosome_names = []
        with open(self.chromosome_name_filename, 'r+') as chromosome_name_file:
            for line in chromosome_name_file:
                chromosome_names.append(ChromosomeName(line.rstrip('\n')))

        # Print the original list, sort and print the sorted list.
        print('\n'.join([chromosome_name.content for chromosome_name in chromosome_names]))
        print('---')
        results = sorted(chromosome_names)
        print('\n'.join([chromosome_name.content for chromosome_name in results]))

    @staticmethod
    def main():
        """
           Entry point into the chromosome sorter.  Parses the argument list, which presently includes just the
           path of the chromosome list file and invokes the sorter.
           """
        parser = argparse.ArgumentParser(description='Sort Chromosome Names')
        parser.add_argument('-c', '--chromosome_name_filename',
                            help="Path to file containing one chromosome name per line of the chromosome"
                                 " to be sorted.")
        args = parser.parse_args()
        print(args)
        chromosome_sort = ChromosomeSort(args.chromosome_name_filename)
        chromosome_sort.sort_chromosome_names()


class ChromosomeName:
    """
    Class holds the chromosome name and provides methods for making comparisons.  Case is disregarded.
    """

    def __init__(self, content):
        """
        Provide a lower case version of the string representation of the chromosome name.
        :param content: string representation of the chromosome name
        """
        self.content = content.lower()

    def __eq__(self, other):
        """
        Test equivalence.  Trivial comparison of the ChromosomeName string with that of another.  None
        ChromosomeName objects are never equivalent.
        :param other: other ChromosomeName object for comparison
        :return: True if equivalent and false otherwise.
        """
        return isinstance(other, ChromosomeName) and self.content == other.content

    def __ne__(self, other):
        """
        Test non-equivalence.  Just the negation of __eq__()
        :param other: other ChromosomeName object for comparison
        :return: True if not equivalent and false otherwise.
        """
        return not self.__eq__(other)

    def __gt__(self, other):
        """
        Where the real comparison work is done.  Makes comparison between this ChomosomeName and the one
        provided based upon the ordering rules described in the ChomosomeSort class documentation.
        :param other: other ChromosomeName object for comparison
        :return: True if this ChromosomeName is greater and false otherwise.
        """

        # Remaining content > no content
        if not self.content:
            return False
        if not other.content:
            return True

        # Names starting with digits > names not starting with digits
        if not self.content[0].isdigit() and other.content[0].isdigit():
            return True
        if self.content[0].isdigit() and not other.content[0].isdigit():
            return False

        # Extract leading digits, if any
        name_match = re.match('^(\d+)(.*)$', self.content)
        other_name_match = re.match('^(\d+)(.*)$', other.content)

        # If both name start with digits order by digits, the name having digits composing the larger
        # number is greater.
        if name_match and other_name_match:
            if int(name_match.group(1)) > int(other_name_match.group(1)):
                return True
            if int(name_match.group(1)) < int(other_name_match.group(1)):
                return False

            # if the digits match, check whether one chromosome name object created from the remaining string is
            # greater than the other chromosome name object.
            return ChromosomeName(name_match.group(2)) > ChromosomeName(other_name_match.group(2))

        # Names starting with roman numerals and special XYM designations > names not starting with roman numerals
        # and special characters
        if not re.match('^([ivxym]+)(.*)$', self.content) and re.match('^([ivxym]+)(.*)$', other.content):
            return True
        if re.match('^([ivxym]+)(.*)$', self.content) and not re.match('^([ivxym]+)(.*)$', other.content):
            return False

        # Extract leading roman numerals, if any
        name_match = re.match('^([ivx]+)(.*)$', self.content)
        other_name_match = re.match('^([ivx]+)(.*)$', other.content)
        if name_match and other_name_match:
            name_roman = Roman(name_match.group(1))
            other_roman = Roman(other_name_match.group(1))
            if name_roman > other_roman:
                return True
            if name_roman < other_roman:
                return False

            # if the roman numerals match, check whether one chromosome name object created from the remaining string
            # is greater than the other chromosome object.
            return ChromosomeName(name_match.group(2)) > ChromosomeName(other_name_match.group(2))

        # X or Y before M
        if self.content == 'm' and (other.content == 'x' or other.content == 'y'):
            return True
        if (self.content == 'x' or self.content == 'y') and other.content == 'm':
            return False

        # Names not starting with 'chr' > names starting with 'chr'
        if not re.match('^chr', self.content) and re.match('^chr', other.content):
            return True
        if re.match('^chr', self.content) and not re.match('^chr', other.content):
            return False

        # Extract leading 'chr' if any
        name_match = re.match('^(chr)(.*)$', self.content)
        other_name_match = re.match('^(chr)(.*)$', other.content)

        # If both chr, check whether one chromosome name object created from the remaining string
        # is greater than the other chromosome object
        if name_match and other_name_match:
            return ChromosomeName(name_match.group(2)) > ChromosomeName(other_name_match.group(2))

        # Names not starting with alphabetic characters > names starting with alphabetic characters.
        if not self.content[0].isalpha() and other.content[0].isalpha():
            return True
        if self.content[0].isalpha() and not other.content[0].isalpha():
            return False

        # Extract leading alphabetic content, if any
        name_match = re.match('^([a-z]+)(.*)$', self.content)
        other_name_match = re.match('^([a-z]+)(.*)$', other.content)
        if name_match and other_name_match:

            # Which character string is greater determine by string >
            if name_match.group(1) > other_name_match.group(1):
                return True
            if name_match.group(1) < other_name_match.group(1):
                return False

            # if both alphabetic strings match, check whether one chromosome name object created from the
            # remaining string is greater than the other chromosome object.
            result = ChromosomeName(name_match.group(2)) > ChromosomeName(other_name_match.group(2))
            return result

        # At this point, the leading characters are not digits, roman numerals, special gender related chars or
        # alphabetic characters, leaving just non-alphanumeric characters.  Since there are no rules for how
        # non-alphanumeric characters are orders we will rely on string comparisons.
        # TODO - add parsing and handling for leading non-alphanumerics
        return True

    def __lt__(self, other):
        """
        Obtained by determining whether the ChromosomeName is greater or equal and negating that finding.
        :param other: other ChromosomeName object for comparison
        :return: True if this ChromosomeName is less and false otherwise.
        """
        return not self.__gt__(other) and self.__ne__(other)


class Roman:
    """
    Class to hold string identified as a roman numeral string, along with its arabic number equivalent.  Comparisons
    are based on a comparison of the associated arabic numbers.
    """

    def __init__(self, content):
        self.content = content
        self.arabic_equivalent = roman.fromRoman(content.upper())

    def __eq__(self, other):
        return isinstance(other, Roman) and self.content == other.content

    def __ne__(self, other):
        return not self.__eq__(other)

    def __gt__(self, other):
        return self.arabic_equivalent > other.arabic_equivalent

    def __lt__(self, other):
        return self.arabic_equivalent < other.arabic_equivalent


if __name__ == "__main__":
    sys.exit(ChromosomeSort.main())
