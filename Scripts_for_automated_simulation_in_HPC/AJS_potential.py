import subprocess


# Parameters to modify
script_path="./job_ShearCoupling_constant_potential.sh"

with open('GB_names.txt') as gn:
    GB_type = [line.rstrip()[:-4] for line in gn][10:11]

#GB_type = ['struct_S7_fcc_N4_1_2_Al_M99_210802.1717','struct_S13a_fcc_N0_n5_1_Al_M99_201228.449','struct_S25a_fcc_N0_1_7_Al_M99_220118.4463','struct_S29b_fcc_N2_7_11_Al_M99_210713.4193']

T_values=[250,350,400,500,550,600]
tau_values=[100,1300,1500,1700,2000,2300]
latP_values = [4.0449,4.0544,4.0594,4.0699,4.0755,4.0812]

def modify_job_script(script_path, parameter1, parameter2, parameter3,parameter4):
    # Read the original job script
    with open(script_path, 'r') as file:
        lines = file.readlines()

    # Modify the parameters in the script
    modified_lines = []
    for line in lines:
        if line.startswith('GB='):
            line = f'GB="{parameter1}"\n'
        elif line.startswith('T='):
            line = f'T={parameter2}\n'
        elif line.startswith('tau='):
            line = f'tau={parameter3}\n'
        elif line.startswith('latP='):
            line = f'latP={parameter4}\n'
        elif line.startswith('#SBATCH -J'):
            line = f'#SBATCH -J {parameter1}_{parameter2}K_{parameter3}Gpa\n'
        elif line.startswith('#SBATCH --output'):
            line = f'#SBATCH --output {parameter1}_{parameter2}K_{parameter3}_potential.out\n'
        modified_lines.append(line)

    # Write the modified script to a temporary file
    modified_script_path = 'modified_ShearCoupling_constant_potential.sh'
    with open(modified_script_path, 'w') as file:
        file.writelines(modified_lines)

    return modified_script_path


def submit_job_script(script_path):
    # Submit the job script to Slurm using sbatch command
    subprocess.run(['sbatch', script_path])



for p1 in GB_type:
    for index,p2 in enumerate(T_values):
        for p3 in tau_values:
            #print(p1)
            modified_script_path = modify_job_script(script_path, p1, p2, p3,latP_values[index])
            submit_job_script(modified_script_path)
