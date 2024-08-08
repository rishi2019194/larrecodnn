import re


def timing_analysis(filename):
    # Initialize a list to store the extracted numbers
    numbers = []

    # Define the regular expression pattern for extracting numbers
    pattern = r'Time taken for inference: (\d+\.?\d*) seconds'

    # Open the file and process it line by line
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                # Extract the number from the match object
                number = float(match.group(1))
                numbers.append(number)
                
    # Calculate the average if there are any numbers
    if numbers:
        average = sum(numbers) / len(numbers)
        print(f'Average time taken for inference: {average:.2f} seconds')
    else:
        print('No numbers found in the file.')

print("Nugraph analysis with CPU & jit model: ")
timing_analysis("testinference_nugraph.txt")
print()

print("EAF-server analysis with GPU(normal triton): ")
timing_analysis("testinference_triton_gpu.txt")
print()

print("EAF-server analysis with CPU(normal triton): ")
timing_analysis("testinference_triton_cpu.txt")
print()

print("EAF-server analysis with GPU(Nusonic triton): ")
timing_analysis("testinference_nusonic_triton_gpu.txt")
print()

print("EAF-server analysis with CPU(Nusonic triton): ")
timing_analysis("testinference_nusonic_triton_cpu.txt")
print()

print("Apptainer-server analysis with CPU(normal triton): ")
timing_analysis("testinference_apptainer_triton.txt")
print()

print("Apptainer-server analysis with CPU(Nusonic triton): ")
timing_analysis("testinference_apptainer_nusonic_triton.txt")
print()
