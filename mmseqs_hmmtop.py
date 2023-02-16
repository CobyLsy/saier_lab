import subprocess
import os


#Ask for names of query and target databases, and Ask to name the result database
query_str = input("Enter query database: ")
target_str = input("Enter target database: ")
result_str = input("Name the result database from mmseqs search: ")

#Perform mmseqs_search
search_cmd = ['mmseqs', 'search', query_str, target_str, result_str, 'tmp']
if os.path.isfile(result_str + ".dbtype") == False:
    subprocess.run(search_cmd, check = True)

#Perform mmseqs_convertalis
flags = 'query,target,alnlen,evalue,bits,pident,qstart,qend,qlen,qcov,tstart,tend,tlen,tcov'
convertalis_cmd = ['mmseqs', 'convertalis', query_str, target_str, result_str, result_str + ".tab", "--format-mode", "4","--format-output", flags]
if os.path.isfile(result_str + ".tab") == False:
    subprocess.run(convertalis_cmd, check = True)

#Perform hmmtop
faa_str = input('Enter file name to run hmmtop on: ')
hmmResult_str = input('Name the result file from hmmtop (name only, without .hmmtop): ')
hmmtop_cmd = ['hmmtop', '-if=' + faa_str, '-of=' + hmmResult_str + '.hmmtop']
if os.path.isfile(hmmResult_str + '.hmmtop') == False:
    subprocess.run(hmmtop_cmd, check = True)

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
        # Create a dictionary using the header as the keys and the fields as the values
        blast_result = {header[i]: fields[i] for i in range(len(header))}
        # Use the first field as the key to the nested dictionary
        key = blast_result.pop(header[0])
        # Add the blast result to the nested dictionary
        mmseqs_dict[key] = blast_result

# Print the resulting dictionary
print(mmseqs_dict)

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
        # ProteinLen indication is parsed into two fields because of the space
        # restore protein length to one object
        combined_protlen = fields[0] + fields[1]
        # create list to store pairs of tms
        tms_list = []
        for i in range(5, len(fields), 2):
            tms_pairs = [fields[i], fields[i+1]]
            tms_list.append(tms_pairs)
        hmmtop_result = {fieldNames[0]:combined_protlen, fieldNames[1]:fields[4], fieldNames[2]:tms_list}
        hmmtop_dict[fields[2]] = hmmtop_result
    

# Print the resulting dictionary
print(hmmtop_dict)
