with open('text.txt','r') as test:
    list_of_lines = [line for line in test.readlines()]
    rmsd=[x[1:3] for x in list_of_lines]
    print(rmsd)



    # 19:24