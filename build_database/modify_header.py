import re
import os

files = os.listdir(os.getcwd())
for f in files:
    if f.endswith('faa'):
        print(f)
        with open(f, 'r') as f1:
            f_modified_path = f.replace('.faa','_modified.faa')
            if os.path.exists(f_modified_path):
                os.remove(f_modified_path)
            with open(f_modified_path, 'a') as f2:
                i = 1
                line = f1.readline()
                while (line):
                    if line[0] == '>':
                        print('That is the header!')
                        f2.write(">gene_%d\n" % i)
                        print(">gene_%d\n" % i)
                        i += 1
                    else:
                        f2.write(line)
                        #print(line)
                    line = f1.readline()