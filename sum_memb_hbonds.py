import os

common_name = 'outfit'
matching_files = [filename for filename in os.listdir() if common_name in filename and filename.endswith(".txt")]
file_last_column_sum = {}
for filename in matching_files:
    input_file_path = os.path.join(os.getcwd(), filename)
    with open(input_file_path, 'r') as input_file:
        last_column_sum = 0.0
        for line in input_file:
            columns = line.strip().split()
            if columns:
                last_column = columns[-1]
                try:
                    value = float(last_column[:-1])
                    last_column_sum += value
                except ValueError:
                    pass 

        file_last_column_sum[filename] = last_column_sum
sorted_results = sorted(file_last_column_sum.items(), key=lambda x: x[1], reverse=True)

#save output
output_file_path = 'sorted_results.txt'
with open(output_file_path, 'w') as output_file:
    for filename, last_column_sum in sorted_results:
        output_file.write(f'{filename}: {last_column_sum:.2f}\n')
print(f"Sorted results saved in the '{output_file_path}' file.")
