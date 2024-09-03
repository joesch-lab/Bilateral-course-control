import numpy as np

def get_expected_response_hemi(opto_response, p=2, translation='N'):
    exp_stim = np.bitwise_or(opto_response[:, 0].astype(np.int), opto_response[:, 1].astype(np.int)) ## assign the expected direction of response, e.g. 1,0 => 1; 0,-1 => -1; 1,-1, => -1
    translation_indices = np.argwhere(np.bitwise_xor(opto_response[:, 0].astype(np.int), opto_response[:, 1].astype(np.int)) == -2).flatten() ## check if the stimulus was translatory, 1, -1 or -1, 1
    left_stronger = np.argwhere(opto_response[translation_indices, 2] > opto_response[translation_indices, 3]).flatten() ## check if contrast or frequency value is higher for left side
    right_stronger = np.argwhere(opto_response[translation_indices, 2] < opto_response[translation_indices, 3]).flatten()
    exp_stim[translation_indices][left_stronger] = opto_response[translation_indices, 0][left_stronger]
    exp_stim[translation_indices][right_stronger] = opto_response[translation_indices, 0][right_stronger]
    print(p)
    opto_response[:, p] = np.array(exp_stim) ## the first entry is the expected direction
    Ftb_indices = np.argwhere(((opto_response[:, 0]==-1) & (opto_response[:, 1]==0)) |
                              ((opto_response[:, 0]==0) & (opto_response[:, 1]==1))).flatten()
    BtF_indices = np.argwhere(((opto_response[:, 0]==1) & (opto_response[:, 1]==0)) |
                              ((opto_response[:, 0]==0) & (opto_response[:, 1]==-1))).flatten()
    full_rotation_indices = np.argwhere(((opto_response[:, 0]==1) & (opto_response[:, 1]==1)) |
                              ((opto_response[:, 0]==-1) & (opto_response[:, 1]==-1))).flatten()

    exp_type = np.ones(exp_stim.shape[0])*5
    exp_type[Ftb_indices] = -1
    exp_type[BtF_indices] = 0
    exp_type[full_rotation_indices] = 1
    if translation == 'Y':
        fwd_indices = np.argwhere(((opto_response[:, 0] == -1) & (opto_response[:, 1] == 1))).flatten()
        bwd_indices = np.argwhere(((opto_response[:, 0] == 1) & (opto_response[:, 1] == -1))).flatten()
        exp_type[fwd_indices] = 2
        exp_type[bwd_indices] = 3
    ## create labels for the 3 conditions
    labels = []
    if translation == 'Y':
        for i in [-1, 0, 1, 2, 3]:
            if i == -1:
                labels.append('Front to Back')
            elif i == 0:
                labels.append('Back to Front')
            elif i == 1:
                labels.append('Full rotation')
            elif i == 2:
                labels.append('FtB translation')
            elif i == 3:
                labels.append('BtF translation')
    else:
        for i in [-1, 0, 1]:
            if i == -1:
                labels.append('Front to Back')
            elif i == 0:
                labels.append('Back to Front')
            elif i == 1:
                labels.append('Full rotation')
    return opto_response, exp_type, labels
