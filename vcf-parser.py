import csv


def clinvar_vcf_to_tsv(clinvar_vcf_file_name):
    """

    :param :
    :return :
    """

    # Loops through the vcf file line by line
    with open(clinvar_vcf_file_name) as vcf_file:
        vcf_list = []
        for line in vcf_file:
            templist = []
            # If the line does not start with an # it appends the
            # variant information to a list
            if not line.startswith("#"):
                templist.append(line.split("\t")[0])  # Chromosome
                templist.append(line.split("\t")[1])  # Position
                templist.append(line.split("\t")[3])  # REF
                templist.append(line.split("\t")[4])  # ALT

                info_col = line.split("\t")[7]

                # Gene name + ID (format name:id,
                # (When it contains multiple genes: name1:id1|name2:id2)
                if "GENEINFO=" in info_col:
                    templist.append(
                        info_col.split("GENEINFO=")[1].split(";")[0])
                else:
                    templist.append("null")

                templist.append("Wild_type")  # Wild type?

                # Variant type
                if "MC=" in info_col:
                    templist.append(
                        info_col.split("MC=")[1].split(";")[0].split("|")[1])
                else:
                    templist.append("null")

                # Clinical significance
                if "CLNSIG=" in info_col:
                    templist.append(info_col.split("CLNSIG=")[1].split(";")[0])
                else:
                    templist.append("null")

                # Review status
                if "CLNREVSTAT=" in info_col:
                    templist.append(
                        info_col.split("CLNREVSTAT=")[1].split(";")[0])
                else:
                    templist.append("null")

                vcf_list.append(templist)

            # Loops through first x amount of  variants
            if len(vcf_list) == 20000:
                break

    vcf_file.close()

    with open("clinvar_output.tsv", 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(["Chromosome", "Position", "REF", "ALT",
                             "Gene_name", "Wild_type?", "Variant_type",
                             "Clinical_significance", "Review_status"])
        for variant in vcf_list:
            tsv_writer.writerow(variant)

    # Closes the file
    out_file.close()


def main():
    clinvar_vcf_file_name = "clinvar.vcf"
    vcf_file_name = "gnomad.genomes.v3.1.2.sites.chrY.vcf"

    clinvar_vcf_to_tsv(clinvar_vcf_file_name)


main()
