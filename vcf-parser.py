import csv
import os.path
from Bio.SeqUtils import seq3


def downsize_clinvar(clinvar_vcf_file_name):
    with open(clinvar_vcf_file_name) as vcf_file, open(
            "clinvar_small.vcf", "w") as out_file:
        for line in vcf_file:
            if os.path.getsize("clinvar_small.vcf") < 51380224:  # 49MB
                try:
                    if line.startswith("#"):
                        out_file.write(line)
                    elif "CLNSIG" in line.split("\t")[7] and (
                            line.split("\t")[7].split("CLNSIG=")[1].split(";")[
                                0] == "Benign" or
                            line.split("\t")[7].split("CLNSIG=")[1].split(";")[
                                0] == "Pathogenic"):
                        out_file.write(line)
                except:
                    pass
    out_file.close()


def clinvar_vcf_to_tsv(clinvar_vep_standaard_param):
    """
    Writes information from the vep annotated clinvar vcf file to a tsv file
    :param clinvar_vep_standaard_param:
    :return clinvar_output.tsv: TSV file with the columns "Chromosome",
        "Position", "REF", "ALT", "Symbol", "Gene", "Protein_Change", "CLNSIG",
        "PolyPhen", "Consequence".
    """

    # Loops through the vcf file line by line
    with open(clinvar_vep_standaard_param) as vcf_file:
        vcf_list = []
        for line in vcf_file:
            # Creates a temporary list for every new line
            templist = ["Chromosome", "Position", "REF", "ALT", "Symbol",
                        "Gene", "Protein_Change", "CLNSIG", "PolyPhen",
                        "Consequence"]
            # If the line does not start with an # it appends the
            # variant information to a list
            if not line.startswith("#"):
                templist[0] = line.split("\t")[0].replace("chr", "")  # Chr
                templist[1] = line.split("\t")[1]  # Position
                templist[2] = line.split("\t")[3]  # REF
                templist[3] = line.split("\t")[4]  # ALT

                # Splits and saves the info column of the vcf file to a
                # variable
                info_col = line.split("\t")[7]

                # Checks if the CLNSIG is in the info column
                if "CLNSIG=" in info_col:
                    templist[7] = info_col.split("CLNSIG=")[1].split(";")[0]
                else:
                    templist[7] = "null"

                # Splits the vep annotated list to single variants
                vep_list = info_col.split("CSQ=")[1].split(",")

                # Loops through the different vep variants and saves it
                # to a list
                for vep in vep_list:
                    # SYMBOL 3, Gene 4, HGVSp 11, PolyPhen 28,
                    # Consequence 1, Protein_change 15, Protein_pos 14
                    templist[4] = vep.split("|")[3]
                    templist[5] = vep.split("|")[4]
                    templist[8] = vep.split("|")[28]
                    templist[9] = vep.split("|")[1]
                    if not "frameshift_variant" in vep.split("|")[1]:
                        if vep.split("|")[14] and vep.split("|")[15] and "/" \
                                in vep.split("|")[15]:
                            templist[6] = seq3(vep.split("|")[15],
                                               undef_code=vep.split("|")[14])
                            vcf_list.append(templist)
                    else:
                        templist[6] = vep.split("|")[15].replace("/",
                                                                 vep.split(
                                                                     "|")[14])
                        vcf_list.append(templist)

    # Closes the vcf file
    vcf_file.close()

    # Creates a tsv file and writes the retrieved info from the vcf file
    # in it
    with open("clinvar_output.tsv", 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(["Chromosome", "Position", "REF", "ALT", "Symbol",
                             "Gene", "Protein_Change", "CLNSIG", "PolyPhen",
                             "Consequence"])
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
            templist = ["Chromosome", "Position", "REF", "ALT", "Symbol",
                        "Gene", "HGVSp", "PolyPhen", "Consequence"]
            # If the line does not start with an # it appends the
            # variant information to a list
            if not line.startswith("#"):
                # Chromosome
                templist[0] = line.split("\t")[0].replace("chr", "")
                templist[1] = line.split("\t")[1]  # Position
                templist[2] = line.split("\t")[3]  # REF
                templist[3] = line.split("\t")[4]  # ALT

                info_col = line.split("\t")[7]

                # vep
                # SYMBOL 3, Gene 4, HGVSp (ENSP00000372547.1:p.Tyr96Asn) 11, PolyPhen 34, Consequence 1
                if "vep=" in info_col:
                    vep_list = info_col.split("vep=")[1].split(",")
                    for vep in vep_list:
                        try:
                            if vep.split("|")[11] and (
                                    "benign" in vep.split("|")[
                                34] or "probably_damaging" in vep.split("|")[
                                        34]):
                                templist[4] = vep.split("|")[3]
                                templist[5] = vep.split("|")[4]
                                templist[6] = vep.split("|")[11]
                                templist[7] = vep.split("|")[34]
                                templist[8] = vep.split("|")[1]
                                vcf_list.append(templist)
                                break
                        except:
                            pass
                # Loops through first x amount of  variants
                if len(vcf_list) >= 100:
                    break

    vcf_file.close()

    with open("gnomad_output.tsv", 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(["Chromosome", "Position", "REF", "ALT",
                             "Symbol", "Gene", "HGVSp", "PolyPhen",
                             "Consequence"])
        for variant in vcf_list:
            tsv_writer.writerow(variant)

    # Closes the file
    out_file.close()


def main():
    gnomad_vcf_file_name = "gnomad.genomes.v3.1.2.sites.chrY.vcf"
    clinvar_vep_standaard_param = "clinvar_vep_standaard_param.vcf"

    # print("Type '1' if you would like to run clinvar_to_tsv. \n"
    #       "Type '2' if you would like to run gnomad_to_tsv.")
    selected_function = input(
        "Type '1' if you would like to run clinvar_to_tsv. \n"
        "Type '2' if you would like to run gnomad_to_tsv.\n"
        "Type '3' if you would like to run both.\n")

    if selected_function == "1":
        clinvar_vcf_to_tsv(clinvar_vep_standaard_param)
        print("clinvar done")

    elif selected_function == "2":
        gnomad_vcf_to_tsv(gnomad_vcf_file_name)
        print("gnomad klaar")

    elif selected_function == "3":
        clinvar_vcf_to_tsv(clinvar_vep_standaard_param)
        print("clinvar done")

        gnomad_vcf_to_tsv(gnomad_vcf_file_name)
        print("gnomad klaar")

    # Klein clinvar bestand aanmaken:
    # clinvar_vcf_file_name = "clinvar.vcf"
    # downsize_clinvar(clinvar_vcf_file_name)
    # print("clinvar vcf downsized")


main()
