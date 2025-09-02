def get_hod_string(hod_params,use_extension=False):
    # hod_keys = ['logM_cut','logM1','sigma','kappa','alpha'] # base hod
    # if use_extension: hod_keys += ['Acent','Bcent','Asat','Bsat','alpha_c','alpha_s','s','s_v','s_p','s_r']

    hod_keys = hod_params.keys()
    out_str = ''
    for i,key in enumerate(hod_keys): 
        if i==0: out_str += '{:s}_{:.2f}'.format(key,hod_params[key])
        else: out_str += '_{:s}_{:.2f}'.format(key,hod_params[key])
    return out_str
