import subprocess
import os
from subprocess import call

#Ask for names of query and target databases, and Ask to name the result database
query_str = input("Enter query database: ")
target_str = input("Enter target database: ")
result_str = input("Name the result database from mmseqs search: ")

#Perform mmseqs_search
search_cmd = ['mmseqs', 'search', query_str, target_str, result_str, 'tmp', '-a']
if os.path.isfile(result_str + ".dbtype") == False:
    subprocess.run(search_cmd, check = True)

#Perform mmseqs_convertalis
flags = 'query,target,alnlen,evalue,bits,pident,qstart,qend,qlen,qcov,tstart,tend,tlen,tcov,qaln,taln'
convertalis_cmd = ['mmseqs', 'convertalis', query_str, target_str, result_str, result_str + ".tab", "--format-mode", "4","--format-output", flags]
if os.path.isfile(result_str + ".tab") == False:
    subprocess.run(convertalis_cmd, check = True)

#Convert Result From mmseqs alignment to fasta format
awkParse = '{print ">"$1; print $15}'
sedOne = '1,2d'
sedTwo = '/^>/! s/-//g'

if os.path.isfile(result_str + '_final.faa') == False:
    call('awk ' + "'" + awkParse + "'" + ' ' + result_str + '.tab > ' + result_str + '.faa', shell = True)
    call('sed ' + "'" + sedOne + "'" + ' ' + result_str + '.faa > ' + result_str + '_f.faa', shell = True)
    call('sed ' + "'" + sedTwo + "'" + ' ' + result_str + '_f.faa > ' + result_str + '_final.faa', shell = True)
    call('rm ' + result_str + '.faa', shell = True)
    call('rm ' + result_str + '_f.faa', shell = True) 

print("Converted result file from mmseqs search to FASTA format as: " + result_str + "_final.faa" + " for use in hmmtop")

#Perform hmmtop
faa_str = input('Enter file name to run hmmtop on: ')
hmmResult_str = input('Name the result file from hmmtop (name only, without .hmmtop): ')
hmmtop_cmd = ['hmmtop', '-if=' + faa_str, '-of=' + hmmResult_str + '.hmmtop']
if os.path.isfile(hmmResult_str + '.hmmtop') == False:
    subprocess.run(hmmtop_cmd, check = True)

#Download Substrate Data
substrate_str = input('Name the file to store substrate data (without .tsv): ')
substrate_cmd = ['wget', '-O', substrate_str + '.tsv', 'https://tcdb.org/cgi-bin/substrates/getSubstrates.py']
if os.path.isfile(substrate_str + '.tsv') == False:
    subprocess.run(substrate_cmd, check = True)

#Open the tab file for reading
with open(result_str + ".tab", 'r') as f:
    # Read the header line and split it by tabs
    header = f.readline().strip().split('\t')
    # Initialize an empty dictionary
    mmseqs_dict = {}
    # Read through each line of the file
    for line in f:
        # Split the line by tabs
        fields = line.strip().split('\t')
        fields[2] = int(fields[2]) #alnlen
        fields[3] = float(fields[3]) #evalue
        fields[4] = int(fields[4]) #bits
        fields[5] = float(fields[5]) #pident
        fields[6] = int(fields[6]) #qstart
        fields[7] = int(fields[7]) #qend
        fields[8] = int(fields[8]) #qlen
        fields[9] = float(fields[9]) #qcov
        fields[10] = int(fields[10]) #tstart
        fields[11] = int(fields[11]) #tend
        fields[12] = int(fields[12]) #tlen
        fields[13] = float(fields[13]) #tcov

        # Create a dictionary using the header as the keys and the fields as the values
        blast_result = {header[i]: fields[i] for i in range(len(header))}
        # Use the first field as the key to the nested dictionary
        key = blast_result.pop(header[0])
        # Add the blast result to the nested dictionary
        mmseqs_dict[key] = blast_result

# Print the resulting dictionary
#print(mmseqs_dict)

# Open the .hmmtop file for reading
with open( hmmResult_str + '.hmmtop', 'r') as f:
    # Initialize an empty dictionary
    hmmtop_dict = {}
    # Manually create names for fields
    fieldNames = ["protlen", "ntms","tms"]
    # Read through each line of the file
    for l in f:
        # Split the line by spaces
        fields = l.strip().split()
        # empty list for individual line
        hmmtop_result = {}
        # create list to store pairs of tms
        tms_list = []
        for i in range(5, len(fields), 2):
            tms_pairs = [int(fields[i]), int(fields[i+1])]
            tms_list.append(tms_pairs)
        hmmtop_result = {fieldNames[0]:int(fields[1]), fieldNames[1]:int(fields[4]), fieldNames[2]:tms_list}
        hmmtop_dict[fields[2]] = hmmtop_result
    

# Print the resulting dictionary
print(hmmtop_dict)

with open( substrate_str + '.tsv', 'r') as f:
    # Initialize an empty dictionary
    substrate_dict = {}
    # Read through each line of the file
    for line in f:
        # Split the line by tabs
        fields = line.strip().split('\t')
        chebi = fields[1].split('|')
        element = []
        for i in chebi:
            rawPair = i.split(';')
            pair = [int(rawPair[0].split(":")[1]), rawPair[1]]
            element.append(pair)
        substrate_dict[fields[0]] = element

#print(substrate_dict)