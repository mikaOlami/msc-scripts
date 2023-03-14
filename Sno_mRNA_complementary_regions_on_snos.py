#!/usr/bin/env python
# Mika Olami 10/01/2023
import sys
import pandas as pd
import xlsxwriter

# create an excel file with charts of the complementary regions of snos-mRNAs on the snoRNAs USAGE: python3
# Sno_mRNA_complementary_regions_on_snos.py [list of snos with complementary regions coords on snos] [snos fasta file
# (uniq)]

#   how to use
if len(sys.argv) < 3:
    exit("\ncreate an excel file with charts of the complementary regions of snos-mRNAs on the snoRNAs\nUSAGE:\n"
         "python3 Sno_mRNA_complementary_regions_on_snos.py [list of snos with complementary regions coords on snos] [snos fasta file (uniq)]\n")

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)

frames = []  # list of all dataframes
names = []   # list of all snoRNAs
snos_regions_dict = {}  # [sno] -> [list of complementary regions]
frames_dict = {}        # [sno] -> dataframe
snos_fasta = {}         # [sno] -> sequence
snos_comp_list_dict = {}  # [sno] -> [list of 0,1,2... #complementary mRNAs in each bp in the sno]

# create dictionary of sequences for each snoRNA
def get_fasta():
    with open(sys.argv[2], "r") as fasta:
        for line in fasta:
            line = line.strip()
            if line.startswith(">"):
                sno = line.split(":")[0].strip(">")
            else:
                snos_fasta[sno] = line


# create dictionary for list of #complementary mRNAs on each bp in the snoRNA
def create_comp_regions():
    for sno in snos_regions_dict:
        snos_comp_list_dict[sno] = []
        length = len(snos_fasta[sno])
        for i in range(length):
            snos_comp_list_dict[sno].append(0)
        for comp_region in snos_regions_dict[sno]:
            start_reg = int(comp_region.split('-')[0])
            end_reg = int(comp_region.split('-')[1])
            for i in range(start_reg - 1, end_reg):
                snos_comp_list_dict[sno][i] += 1


# create dataframe for every snoRNA
def create_df():
    for sno in snos_regions_dict:
        sequence = snos_fasta[sno]
        comp_mRNAs_in_region = snos_comp_list_dict[sno]
        d = {'sequence': [x for x in sequence], '#complementary mRNA partners in region': comp_mRNAs_in_region}
        df = pd.DataFrame(data=d)
        frames.append(df)
        names.append(sno)


# create list of snoRNAs name
def get_names_and_regions():
    with open(sys.argv[1], "r") as input_file:
        for line in input_file:
            line = line.strip()
            l = line.split("\t")
            sno_name = l[0]
            complementary_region = l[1]
            if sno_name not in snos_regions_dict:
                snos_regions_dict[sno_name] = []
            snos_regions_dict[sno_name].append(complementary_region)


# create an excel file of multiple sheets
def create_excel():
    with pd.ExcelWriter('output.xlsx', engine='xlsxwriter') as writer:
        for i in range(len(frames)):
            frames[i].to_excel(writer, sheet_name=names[i], index=False)
            # add chart
            workbook = writer.book
            worksheet = writer.sheets[names[i]]
            chart = workbook.add_chart({'type': 'column'})
            # Get the dimensions of the dataframe.
            (max_row, max_col) = frames[i].shape
            for j in range(1, max_col):
                # Configure the series of the chart from the dataframe data.
                if j == max_col - 1:
                    chart.add_series({'name': [names[i], 0, j], 'categories': [names[i], 1, 0, max_row, 0],
                                      'values': [names[i], 1, j, max_row, j], 'line': {'width': max_row}})
                else:
                    chart.add_series({'name': [names[i], 0, j], 'values': [names[i], 1, j, max_row, j]})
            chart.set_title({'name': names[i] + ' complementary regions to mRNA partners'})
            chart.set_size({'x_scale': 1.5, 'y_scale': 1.25})
            chart.set_y_axis({'name': '#complementary mRNA partners'})
            chart.set_x_axis({'name': 'snoRNA'})
            chart.set_legend({'none': True})
            # Insert the chart into the worksheet.
            worksheet.insert_chart(1, 4, chart)


print("Excel file saved at: output.xlsx")

get_names_and_regions()
get_fasta()
create_comp_regions()
create_df()
create_excel()
