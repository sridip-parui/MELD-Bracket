def get_constraints(fl_name):
     complex = [[] for i in range(len(fl_name))]
     for idx,fl in enumerate(fl_name):
        fl = open("%s" %fl,'r')
        for line in fl:
                #line = line[0:-1].split()
                line = line.strip("\n").split()
                resid1 = int(line[0])
                name1 = line[1]
                resid2 = int(line[2])
                name2 = line[3]
                temp = tuple([resid1,name1,resid2,name2])
                complex[idx].append(temp)
     ##print complex, len(complex)

     return complex
if __name__ == "__main__":
        fl_name=["L00_contacts.txt"]
        print (get_constraints(fl_name))
