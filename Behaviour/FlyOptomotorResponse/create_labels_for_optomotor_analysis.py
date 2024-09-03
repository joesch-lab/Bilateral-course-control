def BTF_FTB_labels(conditions):
    names = []
    if len(conditions[0]) == 1:
        names = [str(x) for x in list(conditions[0][0].astype(int))]
    else:
        for i in range(len(conditions)):
            #left
            if conditions[i][0] == 1:
                new_name = 'BTF'
            elif conditions[i][0] == -1:
                new_name = 'FTB'
            elif conditions[i][0] == 0:
                new_name = 'NO'
            #right
            if conditions[i][1] == 1:
                new_name = new_name + '_' + 'FTB'
            elif conditions[i][1] == -1:
                new_name = new_name + '_' + 'BTF'
            elif conditions[i][1] == 0:
                new_name = new_name + '_' + 'NO'
            names.append(new_name)
    return names


def arrow_labels(conditions):
    names = []
    if len(conditions[0]) == 1:
        names = [str(x) for x in list(conditions[0][0].astype(int))]
    elif len(conditions[0]) == 2:
        for i in range(len(conditions)):
            #left
            if conditions[i][0] == 1:
                new_name = '>'
            elif conditions[i][0] == -1:
                new_name = '<'
            elif conditions[i][0] == 0:
                new_name = 'O'
            else:
                new_name = new_name +' ' + str(conditions[i][0]) + ' \n'
            #right
            if conditions[i][1] == 1:
                new_name = new_name + ' ' + '>'
            elif conditions[i][1] == -1:
                new_name = new_name + ' ' + '<'
            elif conditions[i][1] == 0:
                new_name = new_name + ' ' + 'O'
            else:
                new_name = new_name +' '+str(conditions[i][1])+' \n'
            names.append(new_name)
    elif len(conditions[0]) == 3:
        for i in range(len(conditions)):
            #left
            if conditions[i][0] == 1:
                new_name = '>'
            elif conditions[i][0] == -1:
                new_name = '<'
            elif conditions[i][0] == 0:
                new_name = 'O'
            else:
                new_name = new_name +' '+str(conditions[i][0])+' \n'
            #center
            if conditions[i][2] == 1:
                new_name = new_name + ' ' + '>'
            elif conditions[i][2] == -1:
                new_name = new_name + ' ' + '<'
            elif conditions[i][2] == 0:
                new_name = new_name + ' ' + 'O'
            else:
                new_name = new_name +' '+str(conditions[i][2])+' \n'
            #right
            if conditions[i][1] == 1:
                new_name = new_name + ' ' + '>'
            elif conditions[i][1] == -1:
                new_name = new_name + ' ' + '<'
            elif conditions[i][1] == 0:
                new_name = new_name + ' ' + 'O'
            else:
                new_name = new_name +' '+str(conditions[i][1])+' \n'
            names.append(new_name)
    elif len(conditions[0]) == 4:
        for i in range(len(conditions)):
            #center
            if conditions[i][2] == 1:
                new_name = '  >  \n'
            elif conditions[i][2] == -1:
                new_name = '  <  \n'
            elif conditions[i][2] == 0:
                new_name = '  O  \n'
            else:
                new_name = ' '+str(int(conditions[i][2]))+' \n'
            #left
            if conditions[i][0] == 1:
                new_name = new_name + '>'
            elif conditions[i][0] == -1:
                new_name = new_name + '<'
            elif conditions[i][0] == 0:
                new_name = new_name + 'O'
            else:
                new_name = new_name +' '+str(int(conditions[i][0]))+' \n'
            #right
            if conditions[i][1] == 1:
                new_name = new_name + '   >\n'
            elif conditions[i][1] == -1:
                new_name = new_name + '   <\n'
            elif conditions[i][1] == 0:
                new_name = new_name + '   O\n'
            else:
                new_name = new_name +' '+str(int(conditions[i][1]))+' \n'
            #rear
            if conditions[i][3] == 1:
                new_name = new_name + '  >  '
            elif conditions[i][3] == -1:
                new_name = new_name + '  <  '
            elif conditions[i][3] == 0:
                new_name = new_name + '  O  '
            else:
                new_name = new_name +' '+str(int(conditions[i][3]))+' \n'
            names.append(new_name)
    return names


def CW_CCW_labels(conditions):
    names = []
    if len(conditions[0]) == 1:
        names = [str(x) for x in list(conditions[0][0].astype(int))]
    else:
        for i in range(len(conditions)):
            #left
            if conditions[i][0] == 1:
                new_name = 'CW'
            elif conditions[i][0] == -1:
                new_name = 'CCW'
            elif conditions[i][0] == 0:
                new_name = 'NO'
            #right
            if conditions[i][1] == 1:
                new_name = new_name + '_' + 'CW'
            elif conditions[i][1] == -1:
                new_name = new_name + '_' + 'CCW'
            elif conditions[i][1] == 0:
                new_name = new_name + '_' + 'NO'
            names.append(new_name)
    return names
