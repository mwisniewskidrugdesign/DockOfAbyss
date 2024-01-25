with open('text.txt','r') as test:
    list_of_modes = test.read().split('RMSD\n')[-1].split('\n')[0]
    print(list_of_modes)

    # 19:24