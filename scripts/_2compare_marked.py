#!/usr/local/bin/python3.6


name1 = "germline_total.txt"
name2 = "somatic_high_moderate.txt"

def main():
    onlyIn_name1=name1[:8] +"_only"
    onlyIn_name2=name2[:7] +"_only"
    with open (name1,"r") as f:
        data1=f.read().rstrip()
        data1=data1.split('\n')
    with open (name2,"r") as f:
        data2=f.read().rstrip()
        data2=data2.split("\n")        
    setID1 = find_ID(data1)
    setID2 = find_ID(data2)
    print(len(setID1), len(setID2))
    print(setID2)
    union= setID1 | setID2
    both = setID1 & setID2
    only1 = setID1 - setID2
    only2 = setID2 - setID1
    print("both", len(both),onlyIn_name1, len(only1),onlyIn_name2, len(only2), "union", len(union))
    marked=add_name(data1, data2, only1, only2, both, onlyIn_name1, onlyIn_name2)
    with open ("combined_germline_somatic_marked.txt","w") as f:
        for item in marked:
          f.write(item +"\n")

def find_ID(data):
    set_ID = set()
    for line in data:
        field = line.split('\t')
        if "chr" in field[1]: #chr is in either field[1] or field[2]
            POSID=field[1].strip('\"') +str(field[2])
            set_ID.add(POSID)
        elif "chr" in field[2]:
            POSID=field[2].strip('\"') +str(field[3])
            set_ID.add(POSID)
        else:
            print("chr was not found", line)
    return set_ID
    
def add_name(data1, data2, only1, only2, both, onlyIn_name1, onlyIn_name2):
    marked=[]
    for data in [data1,data2]:
        for line in data:
            field = line.split('\t')
            if "chr" in field[1]: #chr is in either field[1] or field[2]
                POSID=field[1] +str(field[2])
            elif "chr" in field[2]:
                POSID=field[2] +str(field[3])
            else:
                POSID=""
            if POSID in both:
                line = line.rstrip() + "\tBoth"
            elif POSID in only1:
                line = line.rstrip() + "\t" + onlyIn_name1
            elif POSID in only2:
                line = line.rstrip() + "\t" + onlyIn_name2
            else:
                line = line.rstrip() + "\tnot found in both files" 
            marked.append(line)
    marked=set(marked)
    print(len(marked))    
    return marked      
    
main()    

         
