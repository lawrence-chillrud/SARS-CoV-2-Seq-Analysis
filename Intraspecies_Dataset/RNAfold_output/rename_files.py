import os
import glob

dir1 = glob.glob("*.ps")
dir2 = glob.glob("*.mfe")
d = dir1 + dir2
for f in d:
    file_name = f
    new_file_name = ''
    count = 0
    i = 0
    while i < len(file_name):
        char = file_name[i]
        if char != '.' or count != 0:
            new_file_name += char
        if char == '.':
            if count == 0:
                i += 1
            count += 1
        i += 1
    os.rename(file_name, new_file_name)
