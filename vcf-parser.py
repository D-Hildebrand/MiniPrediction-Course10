"""
This module saves the variants in the vcf file to a list.
"""


def vcf_to_list(vcf_file_name):
    """
    Saves the variants in the vcf file to a list.
    :param vcf_file_name: Name of the entered vcf file
    :return vcf_list: List with the structure [CHROM, POS, ID, REF, ALT],
        Might have extra columns added like QUAL, FILTER, INFO depending on the
        entered vcf file.
    """

    # Creates an empty list
    vcf_list = []

    # Loops through the vcf file line by line
    with open(vcf_file_name) as file:
        for line in file:

            # If the line does not start with an # it appends the
            # variant information to a list
            if not line.startswith("#"):

                # Checks if the line starts with chr or M\t, incase
                # a different chromosome column format was used
                if line.startswith("chr"):
                    vcf_list.append(line.replace("chr", "", 1).split("\t"))
                elif line.startswith("M\t"):
                    vcf_list.append(line.replace("M", "MT", 1).split("\t"))
                else:
                    vcf_list.append(line.split("\t"))

    # Closes the file
    file.close()

    # Sorts the list by chromosome
    vcf_list.sort()

    # Returns the list
    return vcf_list
