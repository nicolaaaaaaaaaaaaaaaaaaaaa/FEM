import os
import re
os.chdir('C:\\Users\\nicol\\Desktop\\pythonfiles\\FEM')

# Input file name
input_file_name = 'dkemat.txt'
name_matrix = 'dke'


# Output file name
base_name = os.path.splitext(input_file_name)[0]
output_file_name = base_name + '_out.txt'

with open(input_file_name, 'r') as input_file:
    lines = input_file.readlines()


modified_lines = []

next_line = lines[0].lstrip('#').strip()

for i in range(len(lines)-1):

    line = next_line
    next_line = lines[i+1].lstrip('#').strip()

    if line.endswith('*') and next_line.startswith('*'):
        line = line + '*'
        next_line = next_line.lstrip('*')

    if line.endswith('d') and next_line.startswith('ble'):
        line = line + 'ble('
        next_line = next_line.lstrip('ble(')    
    if line.endswith('db') and next_line.startswith('le'):
        line = line + 'le('
        next_line = next_line.lstrip('le(')  
    if line.endswith('dbl') and next_line.startswith('e'):
        line = line + 'e('
        next_line = next_line.lstrip('e()')   
    if line.endswith('dble'):
        line = line + '('
        next_line = next_line.lstrip('(')      

    if line.endswith('E') and next_line.startswith('E'):
        line = line + 'E'
        next_line = next_line.strip('E')  

    if line.endswith('EE') and next_line.startswith('Y'):
        line = line + 'Y'
        next_line = next_line.strip('Y')   
    if line.endswith('E') and next_line.startswith('EE'):
        line = line + 'EY'
        next_line = next_line.strip('EY')  

    if line.endswith('G') and next_line.startswith('XY'):
        line = line + 'XY'
        next_line = next_line.strip('XY') 
    if line.endswith('GX') and next_line.startswith('Y'):
        line = line + 'Y'
        next_line = next_line.strip('Y') 

    if line.endswith('a') and next_line.startswith('a'):
        line = line + 'a'
        next_line = next_line.strip('a') 

    if line.endswith('b') and next_line.startswith('b'):
        line = line + 'b'
        next_line = next_line.strip('b')  
    
    if line.endswith('n') and next_line.startswith('u'):
        line = line + 'u'
        next_line = next_line.strip('u')   

    # 0.0D0
    pattern1 = re.compile(r'0$')
    pattern2 = re.compile(r'^\.[0-9]D[0-9]')
    if pattern1.search(line) and pattern2.search(next_line):
        digits = next_line[0:3]
        line = line + digits
        next_line = next_line.strip(digits)

    pattern3 = re.compile(r'0\.$')
    pattern4 = re.compile(r'^[0-9]D[0-9]')
    if pattern3.search(line) and pattern4.search(next_line):
        digits = next_line[0:2]
        line = line + digits
        next_line = next_line.strip(digits)

    pattern5 = re.compile(r'0\.[0-9]$')
    pattern6 = re.compile(r'^D[0-9]')
    if pattern5.search(line) and pattern6.search(next_line):
        digits = next_line[0:1]
        line = line + digits
        next_line = next_line.strip(digits)

    pattern7 = re.compile(r'0\.[0-9]D$')
    pattern8 = re.compile(r'^[0-9]')
    if pattern7.search(line) and pattern8.search(next_line):
        digits = next_line[0]
        line = line + digits
        next_line = next_line.strip(digits)

    # 0.00D0
    pattern9 = re.compile(r'0$')
    pattern10 = re.compile(r'^\.[0-9][0-9]D[0-9]')
    if pattern9.search(line) and pattern10.search(next_line):
        digits = next_line[0:4]
        line = line + digits
        next_line = next_line.strip(digits)

    pattern11 = re.compile(r'0\.$')
    pattern12 = re.compile(r'^[0-9][0-9]D[0-9]')
    if pattern11.search(line) and pattern12.search(next_line):
        digits = next_line[0:3]
        line = line + digits
        next_line = next_line.strip(digits)

    pattern13 = re.compile(r'0\.[0-9]$')
    pattern14 = re.compile(r'^[0-9]D[0-9]')
    if pattern13.search(line) and pattern14.search(next_line):
        digits = next_line[0:2]
        line = line + digits
        next_line = next_line.strip(digits)

    pattern15 = re.compile(r'0\.[0-9][0-9]$')
    pattern16 = re.compile(r'^D[0-9]')
    if pattern15.search(line) and pattern16.search(next_line):
        digits = next_line[0:1]
        line = line + digits
        next_line = next_line.strip(digits)

    pattern17 = re.compile(r'0\.[0-9][0-9]D$')
    pattern18 = re.compile(r'^[0-9]')
    if pattern17.search(line) and pattern18.search(next_line):
        digits = next_line[0]
        line = line + digits
        next_line = next_line.strip(digits)


    if next_line.startswith(name_matrix) or next_line.isspace() or not line:
        modified_lines.append(line + '\n')
    else:
        modified_lines.append(line + '&\n')


# Write modified lines to the output file
with open(output_file_name, 'w') as output_file:
    output_file.writelines(modified_lines)

print(f'written {name_matrix} matrix in output called {base_name}_output.txt')
