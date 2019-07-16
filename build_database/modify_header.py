import re
import os

def main():
    files = os.listdir(os.getcwd())
    for f in files:
        i = 1
        if f.endswith('_modified.faa') or f.endswith('_modified.fnn'):
            os.remove(f)
        elif f.endswith('.faa'):
            print(f)
            with open(f, 'r') as f1:
                f_modified_path = f.replace('.faa','_modified.faa')
                if os.path.exists(f_modified_path):
                    os.remove(f_modified_path)
                with open(f_modified_path, 'a') as f2:
                    line = f1.readline()  #read line by line
                    while (line):
                        if line[0] == '>':
                            #print('That is the header!')
                            f2.write(">gene_%d\n" % i)
                            print(">gene_%d\n" % i)
                            i += 1
                        else:
                            f2.write(line)
                        line = f1.readline()
        elif f.endswith('.fnn'):
            print(f)
            with open(f, 'r') as f1:
                f_modified_path = f.replace('.fnn','_modified.fnn')
                if os.path.exists(f_modified_path):
                    os.remove(f_modified_path)
                with open(f_modified_path, 'a') as f2:
                    line = f1.readline()  #read line by line
                    while (line):
                        if line[0] == '>':
                            #print('That is the header!')
                            f2.write(">gene_%d\n" % i)
                            print(">gene_%d\n" % i)
                            i += 1
                        else:
                            f2.write(line)
                        line = f1.readline()

if __name__ == '__main__':
    main()