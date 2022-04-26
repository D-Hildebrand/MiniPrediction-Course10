import csv


def clinvar_vcf_to_tsv(clinvar_vcf_file_name):
    """

    :param :
    :return :
    """
    # TODO: Clinvar alleen pathogeen/benign opslaan naar vcf file voor VEP

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


def gnomad_vcf_to_tsv(gnomad_vcf_file_name):
    """

    :param :
    :return :
    """

    # Loops through the vcf file line by line
    with open(gnomad_vcf_file_name) as vcf_file:
        vcf_list = []
        for line in vcf_file:
            templist = []
            # If the line does not start with an # it appends the
            # variant information to a list
            if not line.startswith("#"):
                # Chromosome
                templist.append(line.split("\t")[0].replace("chr", ""))
                templist.append(line.split("\t")[1])  # Position
                templist.append(line.split("\t")[3])  # REF
                templist.append(line.split("\t")[4])  # ALT

                info_col = line.split("\t")[7]


                #TODO: vep polyphen filteren alleen pathogeen/benign

                # vep
                # SYMBOL 3, Gene 4, VARIANT_CLASS 20, HGVSp (ENSP00000372547.1:p.Tyr96Asn) 11, PolyPhen 34, Consequence 1
                if "vep=" in info_col:
                    vep_list = info_col.split("vep=")[1].split(",")
                    for vep in vep_list:
                        try:
                            if vep.split("|")[11] and vep.split("|")[34]:
                                if len(templist) == 4:
                                    templist.append(vep.split("|")[3])
                                    templist.append(vep.split("|")[4])
                                    templist.append(vep.split("|")[20])
                                    templist.append(vep.split("|")[11])
                                    templist.append(vep.split("|")[34])
                                    templist.append(vep.split("|")[1])
                                else:
                                    templist[4] = vep.split("|")[3]
                                    templist[5] = vep.split("|")[4]
                                    templist[6] = vep.split("|")[20]
                                    templist[7] = vep.split("|")[11]
                                    templist[8] = vep.split("|")[34]
                                    templist[9] = vep.split("|")[1]
                                vcf_list.append(templist)
                        except:
                            pass
                # Loops through first x amount of  variants
                if len(vcf_list) > 100:
                    quit()

    vcf_file.close()

    with open("gnomad_output.tsv", 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(["Chromosome", "Position", "REF", "ALT",
                             "Symbol", "Gene", "Variant class", "HGVSp",
                             "PolyPhen", "Consequence"])
        for variant in vcf_list:
            tsv_writer.writerow(variant)

    # Closes the file
    out_file.close()


def main():
    clinvar_vcf_file_name = "clinvar.vcf"
    gnomad_vcf_file_name = "gnomad.genomes.v3.1.2.sites.chrY.vcf"

    # clinvar_vcf_to_tsv(clinvar_vcf_file_name)
    gnomad_vcf_to_tsv(gnomad_vcf_file_name)


main()
