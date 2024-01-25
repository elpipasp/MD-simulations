import os
import re


file_list = os.listdir(".")

def custom_sort_key(file_name):
    match = re.match(r'saltbr-[A-Z]{3}(\d+)_chain[A-Z]_segname[A-Z]+-[A-Z]{3}(\d+)_chain[A-Z]_segname[A-Z]+\.dat', file_name)
    if match:
        first_number = int(match.group(1))
        second_number = int(match.group(2))
        return first_number, second_number
    return float('inf'), float('inf')
filtered_files = [file_name for file_name in file_list if re.match(r'saltbr-[A-Z]{3}\d+_chain[A-Z]_segname[A-Z]+-[A-Z]{3}\d+_chain[A-Z]_segname[A-Z]+\.dat', file_name)]
sorted_files = sorted(filtered_files, key=custom_sort_key)
for file_name in sorted_files:
    print(file_name)
