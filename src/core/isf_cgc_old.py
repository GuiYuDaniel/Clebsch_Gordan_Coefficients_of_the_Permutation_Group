#
#
#
#
#
#
#
#
# cal ISF and CG of S3~S7
#
#
import os,pickle,copy
import numpy as np

# input
n = 6
r_path = '/Users/guiyu/Desktop/CG/r/'
order_path = '/Users/guiyu/Desktop/CG/CG_order/'
lamada_path = '/Users/guiyu/Desktop/CG/eigenvalue/'
branchlow_path = '/Users/guiyu/Desktop/CG/BranchLaw/'
CG_up_path = '/Users/guiyu/Desktop/CG/CG/'+'S'+str(n-1)+'/'
CG_path = '/Users/guiyu/Desktop/CG/CG/'+'S'+str(n)+'/'
work_path = '/Users/guiyu/Desktop/CG/ISF/'+'S'+str(n)+'/'
Mm_path = '/Users/guiyu/Desktop/CG/YangPic/'
in_path = '/Users/guiyu/Desktop/CG/YamanouchiMatrix(in)/'
lamada_get_r_path = '/Users/guiyu/Desktop/CG/eigenvalue2r/'

# char
n_up = n-1
# functions
# 1
def get_r(n):  # 得到Sn所有构型
    My_file_r = open(r_path + 'S' + str(n) + '.pkl', 'rb')
    r_matrix = pickle.load(My_file_r)
    My_file_r.close()
    return r_matrix

def get_order_r(order_path,n,sigema,miu,r_matrix): # 得到CG序列。方法：以r_matrix去遍历，留下不为零的order
    result_r = []
    result_beta = []
    for r in r_matrix:
        My_file_order = open(order_path + 'S' + str(n) + '/' + str(sigema)+'_'+str(miu)+'/'+str(r)+'.pkl', 'rb')
        order = pickle.load(My_file_order)
        My_file_order.close()
        if order != 0:
            result_r.append(r)
            result_beta.append(order)
    return result_r

def get_lamada_r(n,r):
    path = lamada_path+'/S'+str(n)+'/'+str(r)+'.pkl'
    My_file = open(path,'rb')
    _lamada = pickle.load(My_file)
    My_file.close()
    return _lamada

def get_BranchLow_r(branchlow_path,n,r): # 得到Sn，构型r的分支律
    My_file_BL = open(branchlow_path + 'S' + str(n) +'/BeforeS/'+str(r)+ '.pkl', 'rb')
    r_matrix = pickle.load(My_file_BL)
    My_file_BL.close()
    return r_matrix

def from_lamada_rup_get_r(_lamada,_r_up,n):
    path = lamada_get_r_path+'S'+str(n)+'/'+str(_lamada)+'.pkl'
    My_file = open(path, 'rb')
    _r_matrix = pickle.load(My_file)
    #print('_r_matrix = pickle.load(My_file)',_r_matrix)
    My_file.close()
    if len(_r_matrix) == 1:
        _r = _r_matrix[0]
        return _r
    else:
        for _r_candidate in _r_matrix:
            #print('branchlow_path,n,_r_candidate=',branchlow_path,n,_r_candidate)
            _r_up_candidate_matrix = get_BranchLow_r(branchlow_path,n,_r_candidate)
            if _r_up in _r_up_candidate_matrix:
                return _r_candidate

def get_CG(path):
    My_file_CG = open(path,'rb')
    CG_dict = pickle.load(My_file_CG)
    My_file_CG.close()
    return CG_dict

def get_Mm(path,r,n):
    My_file = open(path+'S'+str(n)+'/'+str(r)+'/'+'Mm.txt','r')
    _m = My_file.readline()
    My_file.close()
    return int(_m)

def calc_fai_matrix(r_sigema_up,r_miu_up,rou_up): # 新的fai加个beta小尾巴
    fai_matrix = []
    for sigema_up in r_sigema_up:
        for miu_up in r_miu_up:
            if n == 3:
                path = CG_up_path + 'sigema' + str(sigema_up) + '_miu' + str(miu_up) + '/rou' + str(
                    rou_up) + '_beta' + str(1) + '_m' + str(1) + '.pkl'
                if os.path.exists(path):  # 搜到了由σ'*μ'得到的结果中有ν'
                    fai_matrix.append([sigema_up, miu_up,1])
            else:
                My_file_order = open(
                    order_path + 'S' + str(n_up) + '/' + str(sigema_up) + '_' + str(miu_up) + '/' + str(rou_up) + '.pkl', 'rb')
                order = pickle.load(My_file_order)
                My_file_order.close()
                if order != 0:
                    for oo in range(1,order+1):
                        fai_matrix.append([sigema_up, miu_up, oo])
    print('fai_matrix=',fai_matrix)
    return fai_matrix

def from_r_rup_mup_get_m(r,r_up,m_up,n=n,n_up=n_up): # 例如由σ,σ',m'得到m
    m = 0
    _n = n
    _n_up = n_up
    bl_r_matrix = get_BranchLow_r(branchlow_path,_n,r)
    for bl_r in bl_r_matrix:
        if r_up != bl_r:
            m_part = get_Mm(Mm_path,bl_r,_n_up)
            m += m_part
        else:
            m += m_up
            return m

def from_in_get_matrix_element(i,r,m1,m2):
    opera_in = [i,n]
    path = in_path +'S'+str(n)+'/in'+str(opera_in)+'/'+str(r)+'.pkl'
    My_file = open(path,'rb')
    in_matrix = pickle.load(My_file)
    My_file.close()
    in_matrix_element = in_matrix[m1-1][m2-1]
    return in_matrix_element

def get_fai_info(fai,sigema,miu,rou_up,beta_up):
    # fai_info = {():int,
    #             ():int,
    #             ...,
    #             'N':int}
    fai_info ={}
    sigema_up = fai[0]
    miu_up = fai[1]
    path = CG_up_path+'sigema'+str(sigema_up)+'_miu'+str(miu_up)+'/rou'+str(rou_up)+'_beta'+str(beta_up)+'_m'+str(1)+'.pkl'
    CG_dict = get_CG(path)
    fai_dict = copy.deepcopy(CG_dict)
    fai_dict.pop('N')
    fai_info['N'] = CG_dict['N']
    #fai_info['lenth'] = len(fai_dict)
    for m1_up,m2_up in fai_dict: # fai_dict.key是（m1，m2）形式，直接call出了m1，m2
        m1 = from_r_rup_mup_get_m(sigema,sigema_up,m1_up)
        m2 = from_r_rup_mup_get_m(miu,miu_up,m2_up)
        fai_info[(m1,m2)] = fai_dict[(m1_up,m2_up)]
    return fai_info

def calc_ISF_matrix_element(fai_1,fai_2,sigema,miu,rou_up,beta_up_1,beta_up_2): # fai是[|σ,σ'>|μ,μ'>]，σ,μ在循环中读入，fai记录的是[σ',μ']
    matrix_element = 0
    for i in range(1,n):
        part_matrix_element = 0
        fai1_info = get_fai_info(fai_1,sigema,miu,rou_up,beta_up_1)
        fai2_info = get_fai_info(fai_2,sigema,miu,rou_up,beta_up_2)
        fai1_N = fai1_info['N']
        fai2_N = fai2_info['N']
        fai1_info.pop('N')
        fai2_info.pop('N')
        for sigema_up_m1,miu_up_m1 in fai1_info:
            for sigema_up_m2,miu_up_m2 in fai2_info:
                #CG_left = np.sqrt(1/fai1_N) * np.sign(fai1_info[(sigema_m1,miu_m1)]) * np.sqrt(fai1_info[(sigema_m1,miu_m1)])
                #CG_right = np.sqrt(1/fai2_N) * np.sign(fai2_info[(sigema_m2,miu_m2)]) * np.sqrt(fai2_info[(sigema_m2,miu_m2)])
                #CG_left = np.sqrt(1/fai1_N)*np.sign(fai1_info[(sigema_up_m1, miu_up_m1)])*np.sqrt(abs(fai1_info[(sigema_up_m1, miu_up_m1)]))
                #CG_right = np.sqrt(1/fai2_N)*np.sign(fai2_info[(sigema_up_m2, miu_up_m2)])*np.sqrt(abs(fai2_info[(sigema_up_m2, miu_up_m2)]))
                CG_left_right = np.sqrt(1/(fai1_N*fai2_N))*np.sqrt(abs(fai1_info[(sigema_up_m1, miu_up_m1)]*fai2_info[(sigema_up_m2, miu_up_m2)]))
                CG_left_right = CG_left_right * np.sign(fai1_info[(sigema_up_m1, miu_up_m1)]*fai2_info[(sigema_up_m2, miu_up_m2)])
                sigema_element = from_in_get_matrix_element(i,sigema,sigema_up_m1,sigema_up_m2)
                miu_element = from_in_get_matrix_element(i,miu,miu_up_m1,miu_up_m2)
                #print('fai1,fai2,四个乘法因子',fai1_info,fai2_info,CG_left,CG_right,sigema_element,miu_element,"in's i=",i)
                #part_matrix_element += CG_left*CG_right*sigema_element*miu_element
                part_matrix_element += CG_left_right * sigema_element * miu_element
        matrix_element += part_matrix_element
    return matrix_element

# TODO 优化解析fai，未来要加上beta操作，未来可能要调整CG系数相位的问题
def calc_ISF_matrix(fai_matrix,sigema,miu,rou_up): # 这里面会因为循环多次解析重复的fai，没关系，先这么做，不行再优化
    _l = len(fai_matrix)
    print('_l=',_l)
    _ISF_matrix = np.zeros([_l,_l])
    for i in range(_l): # 主对角元
        fai_matrix_1 = [fai_matrix[i][0],fai_matrix[i][1]]
        beta_up_1 = fai_matrix[i][2]
        _ISF_matrix[i,i] = calc_ISF_matrix_element(fai_matrix_1,fai_matrix_1,sigema,miu,rou_up,beta_up_1,beta_up_1)
    if _l >= 2:
        for i in range(_l):
            for j in range(i+1,_l):
                fai_matrix_1 = [fai_matrix[i][0], fai_matrix[i][1]]
                beta_up_1 = fai_matrix[i][2]
                fai_matrix_2 = [fai_matrix[j][0], fai_matrix[j][1]]
                beta_up_2 = fai_matrix[j][2]
                _ISF_matrix[i,j] = calc_ISF_matrix_element(fai_matrix_1,fai_matrix_2,sigema,miu,rou_up,beta_up_1,beta_up_2)
                _ISF_matrix[j,i] = _ISF_matrix[i,j]
    print('_ISF_matrix=',_ISF_matrix)
    return _ISF_matrix

def cut_faimatrix_beta(fai_matrix):
    cut_fai_matrix = []
    for ff in fai_matrix:
        cc = [ff[0],ff[1]]
        if not cc in cut_fai_matrix:
            cut_fai_matrix.append([ff[0],ff[1]])
    return cut_fai_matrix

def calc_betas_max_dict(eigenvalues_int):
    '''
    _dict = {lamada:[number,[indexs]]}
    _counter = {lamada:0}
    '''
    _dict = {}
    _counter = {}
    for _index in range(len(eigenvalues_int)):
        _lamada = eigenvalues_int[_index]
        if _lamada in _dict:
            _dict[_lamada][0] += 1
            _dict[_lamada][1].append(_index)
        else:
            _dict[_lamada] = [1,[_index]]
            _counter[_lamada] = 0
    return _dict,_counter

def sym_ant_ly_and_renew_counter(eigenvalue_int,eigenvectors,lamada_beta_counter_dict,same_lamada_dict):
    '''return and action'''
    sym_or_ant = lamada_beta_counter_dict[eigenvalue_int]
    # action
    lamada_beta_counter_dict[eigenvalue_int] += 1

    indexs = same_lamada_dict[eigenvalue_int][1]
    if sym_or_ant == 0:
        eigenvector = eigenvectors[indexs[0]] + eigenvectors[indexs[1]]
    elif sym_or_ant == 1:
        eigenvector = eigenvectors[indexs[0]] - eigenvectors[indexs[1]]
    eigenvector_guiyi = eigenvector * (1/np.sqrt(sum(eigenvector*eigenvector)))
    return eigenvector_guiyi

def schmitt_ly_and_renew_counter(eigenvalue_int,eigenvectors,lamada_beta_counter_dict,same_lamada_dict,finish_schmitt_dict):
    '''return and action'''
    schmit_number = lamada_beta_counter_dict[eigenvalue_int]
    # action1
    lamada_beta_counter_dict[eigenvalue_int] += 1

    indexs = same_lamada_dict[eigenvalue_int][1]
    if schmit_number == 0: # 第一个全加
        eigenvector = 0
        for ss in range(len(indexs)):
            eigenvector += eigenvectors[indexs[ss]]
    else: # 后面的
        # 例如schmit_number=1，指向indexs中第二个也就是indexs[1]，原vector=eigenvectors[indexs[1]]，前面已经造好的有1个，
        # 存在finish_schmitt_dict[0]中
        # 例如schmit_number=3，指向indexs中第二个也就是indexs[3]，原vector=eigenvectors[indexs[3]]，前面已经造好的有3个，
        # 存在finish_schmitt_dict[0]，finish_schmitt_dict[1]，finish_schmitt_dict[2]中
        eigenvector = eigenvectors[indexs[schmit_number]]
        for ss in range(schmit_number):
            cosin = sum(eigenvector * finish_schmitt_dict[ss])
            eigenvector -= cosin * finish_schmitt_dict[ss]
    eigenvector_guiyi = eigenvector * (1/np.sqrt(sum(eigenvector * eigenvector)))
    # action2
    finish_schmitt_dict[schmit_number] = eigenvector_guiyi  # 记录已经做好eigenvector供自己调用{0：vector1,...}
    return eigenvector_guiyi

def jugde_No_0_branchlaw(rou_new,rou_up):
    r_up_matrix = get_BranchLow_r(branchlow_path, n, rou_new)
    if rou_up == r_up_matrix[0]:
        return 0
    else:
        return 1

def jugde_absolute_phase(eigenvector):
    for vector_element in eigenvector:
        if abs(vector_element) < 0.000001:
            pass
        else:
            return np.sign(vector_element)

def get_first_branchlaw(rou_new):
    r_up_matrix = get_BranchLow_r(branchlow_path, n, rou_new)
    return r_up_matrix[0]

def get_right_ISF_first_no_0(sigema,miu,rou_new,beta,rou_up_right):
    My_file = open(work_path+'sigema'+str(sigema)+'_miu'+str(miu)+'/rou'+str(rou_new)+'_beta'+str(beta)+"/rou'"+str(rou_up_right)+"_beta'1.pkl",'rb')
    right_ISF_dict = pickle.load(My_file)
    My_file.close()
    for right_ISF_key in right_ISF_dict:
        if abs(right_ISF_dict[right_ISF_key]) < 0.0001: # 任意一个都可以去比较，所以挑个大一点的安全
            pass
        else:
            return right_ISF_key,right_ISF_dict[right_ISF_key]

def get_m_now(rou_new,rou_up,m_up_now=1):
    r_up_matrix = get_BranchLow_r(branchlow_path, n, rou_new)
    m_now = m_up_now
    for r_son in r_up_matrix:
        if rou_up == r_son:
            return m_now
        else:
            m_r_son = get_Mm(Mm_path,r_son,n_up)
            m_now += m_r_son

def get_m_right(rou_new,rou_up,rou_up_right,m_up_now=1,m_up_up=1): # 相对相位m'选1之后，两重分支律保证必落在绝对相位里
    r_up_up_matrix = get_BranchLow_r(branchlow_path, n-1, rou_up)
    rou_up_up = r_up_up_matrix[0]
    m_up_right = from_r_rup_mup_get_m(rou_up_right,rou_up_up,m_up_up,n=n-1,n_up=n_up-1)
    m_right = from_r_rup_mup_get_m(rou_new,rou_up_right,m_up_right,n=n,n_up=n_up)
    return m_right,m_up_right

def from_sth_get_CG_up_dict(fai,rou,m):
    _sigema = list(fai[0])
    _miu = list(fai[1])
    _rou = rou
    _beta = fai[2]
    _path = CG_up_path+'sigema'+str(_sigema)+'_miu'+str(_miu)+'/rou'+str(_rou)+'_beta'+str(_beta)+'_m'+str(m)+'.pkl'
    _CG_up_dict = get_CG(_path)
    return _CG_up_dict

def calc_ISF_now(eigenvector,sigema,miu,rou_new,beta,fai_right,rou_up_right,fai_matrix,rou_up,m_right,m_up_right,m_now,m_up_now=1):
    dd = from_in_get_matrix_element(n-1, rou_new, m_now, m_right) # i默认就是n-1
    CM = 0
    #print('rou_up_right=',rou_up_right)
    CG_all_dict_right = from_sth_get_CG_up_dict(fai_right,rou_up_right,m_up_right)
    CG_dict_right = copy.deepcopy(CG_all_dict_right)
    CG_dict_right.pop('N')
    CG_dict_of_N_right = CG_all_dict_right['N']
    print('CG_dict_right=',CG_dict_right)
    m_sigema_miu_up_tuple_right = CG_dict_right.keys()
    for index_now in range(len(fai_matrix)):
        fai_now = fai_matrix[index_now]
        MMM = 0
        ISF_eigenvector_element_now = eigenvector[index_now]
        CG_all_dict_now = from_sth_get_CG_up_dict(fai_now,rou_up,m_up_now) # 这个1是要求m_up_now必为1
        CG_dict_now = copy.deepcopy(CG_all_dict_now)
        CG_dict_now.pop('N')
        CG_dict_of_N_now = CG_all_dict_now['N']
        m_sigema_miu_up_tuple_now = CG_dict_now.keys()
        for sigema_up_m_right,miu_up_m_right in m_sigema_miu_up_tuple_right:
            sigema_m_right = from_r_rup_mup_get_m(sigema,list(fai_right[0]),sigema_up_m_right)
            miu_m_right = from_r_rup_mup_get_m(miu,list(fai_right[1]),miu_up_m_right)
            CG_element_right = CG_dict_right[(sigema_up_m_right,miu_up_m_right)]
            for sigema_up_m_now,miu_up_m_now in m_sigema_miu_up_tuple_now:
                sigema_m_now = from_r_rup_mup_get_m(sigema,list(fai_now[0]),sigema_up_m_now)
                miu_m_now = from_r_rup_mup_get_m(miu,list(fai_now[1]),miu_up_m_now)
                CG_element_now = CG_dict_now[(sigema_up_m_now, miu_up_m_now)]
                d_sigema = from_in_get_matrix_element(n-1, sigema, sigema_m_right, sigema_m_now)
                d_miu = from_in_get_matrix_element(n-1, miu, miu_m_right, miu_m_now)
                MMM += np.sign(CG_element_right * CG_element_now) * np.sqrt(abs(CG_element_right * CG_element_now)) * d_sigema * d_miu * (1/np.sqrt(CG_dict_of_N_right*CG_dict_of_N_now))
        CM += MMM * ISF_eigenvector_element_now
    result = CM *(1/dd)
    return result

def from_r_n_get_beat_max(r,n,sigema,miu):
    if n in [2,3]:
        _beta_max = 1
    else:
        My_file_order = open(order_path+'S'+str(n)+'/'+str(sigema)+'_'+str(miu)+'/'+str(r)+'.pkl','rb')
        _beta_max = pickle.load(My_file_order)
        My_file_order.close()
    return _beta_max

def from_betadict_rounew_calc_beta(beta_dict,lamada_r): # 其实这个地方是可以用lamada的
    if lamada_r in beta_dict:
        _beta = beta_dict[lamada_r] + 1
        beta_dict[lamada_r] = _beta
    else:
        _beta = 1
        beta_dict[lamada_r] = _beta
    return _beta

def ask_CG_get_part_CG(path):
    if os.path.exists(path):
        My_file = open(path,'rb')
        CG_part = pickle.load(My_file)
        My_file.close()
    else:
        CG_part = {'N':1}
    return CG_part

def re_new_CG_part(m_sigema, m_miu,CG_preparation,CG_float):
    '''action'''
    if (m_sigema,m_miu) in CG_preparation:
        CG_preparation[(m_sigema, m_miu)] = CG_preparation[(m_sigema,m_miu)] + CG_float
    else:
        CG_preparation[(m_sigema, m_miu)] = CG_float

def save_ISF_by_path(ISF_file_path,save_dict): # 老规矩，存pickle和txt
    Myfile_ISF = open(ISF_file_path+'.pkl', 'wb')
    pickle.dump(save_dict, Myfile_ISF)
    Myfile_ISF.close()
    Myfile_ISF = open(ISF_file_path+'.txt', 'w')
    Myfile_ISF.write(str(save_dict))
    Myfile_ISF.close()
    print(ISF_file_path+'的ISF系数保存成功！系数是：'+str(save_dict))

def save_CG_by_path(CG_file_path,save_dict): # 老规矩，存pickle和txt
    Myfile_CG = open(CG_file_path+'.pkl', 'wb')
    pickle.dump(save_dict, Myfile_CG)
    Myfile_CG.close()
    Myfile_CG = open(CG_file_path+'.txt', 'w')
    Myfile_CG.write(str(save_dict))
    Myfile_CG.close()
    print(CG_file_path+'的CG系数保存成功！系数是：'+str(save_dict))

# input
'''
n = 3
r_path = '/Users/guiyu/Desktop/CG/r/'
order_path = '/Users/guiyu/Desktop/CG/CG_order/'
lamada_path = '/Users/guiyu/Desktop/CG/eigenvalue/'
branchlow_path = '/Users/guiyu/Desktop/CG/BranchLaw/'
CG_up_path = '/Users/guiyu/Desktop/CG/CG/'+'S'+str(n-1)+'/'
CG_path = '/Users/guiyu/Desktop/CG/CG/'+'S'+str(n)+'/'
work_path = '/Users/guiyu/Desktop/CG/ISF/'+'S'+str(n)+'/'
Mm_path = '/Users/guiyu/Desktop/CG/YangPic/'
in_path = '/Users/guiyu/Desktop/CG/YamanouchiMatrix(in)/'
lamada_get_r_path = '/Users/guiyu/Desktop/CG/eigenvalue2r/'

# char
n_up = n-1
'''
if not os.path.exists(work_path):
    os.makedirs(work_path)
    print('建立目录', work_path)

# main
r_sigema = get_r(n) # 定[σ]
r_miu = r_sigema           # 定[μ]
r_matrix = r_sigema
for sigema in r_sigema:
    for miu in r_miu:
        _folder_path = work_path + 'sigema'+str(sigema) + '_' + 'miu'+str(miu) + '/'
        _CG_folder_path = CG_path+'sigema'+str(sigema) + '_' + 'miu'+str(miu) + '/'
        #_CG_folder_path = CG_path  + str(sigema) + '_'  + str(miu) + '/'
        if os.path.exists(_folder_path):  # 有就跳过
            print('跳过:sigema,miu=',sigema,miu)
        else:  # 没有就算
            os.makedirs(_folder_path)
            os.makedirs(_CG_folder_path)
            print('=' * 30, '建立目录', _folder_path, '=' * 30)
        r_rou_up = get_r(n_up)
        r_rou = get_order_r(order_path, n, sigema, miu, r_matrix)  # 定对于σ，μ所有可能的[ν]
        for rou_up in r_rou_up: # 反写循环是为了配合绝对相位约定
            for rou in r_rou:
                if rou_up in get_BranchLow_r(branchlow_path,n,rou):
                    folder_path = _folder_path + 'rou' + str(rou) + '_beta' + str(1) + '/' + "rou'" + str(
                        rou_up) + "_beta'" + str(1) + '.pkl'
                    if os.path.exists(folder_path):  # 有就跳过
                        print("跳过:sigema,miu,rou,rou'=",sigema,miu,rou,rou_up)
                    else:  # 没有就算
                        # TODO 此处一组[ν]——[ν']不仅能算一组[ν]——[ν']哦，可以把[x]——[ν']都算了，所以建立文件夹应该往后放
                        # 下面导入[σ'][μ']
                        r_sigema_up = get_BranchLow_r(branchlow_path, n, sigema)
                        r_miu_up = get_BranchLow_r(branchlow_path, n, miu)
                        fai_matrix = calc_fai_matrix(r_sigema_up, r_miu_up, rou_up)  # [[[σ'],[μ'],β],...]
                        beta_up_max = from_r_n_get_beat_max(rou_up, n_up, fai_matrix[0][0], fai_matrix[0][1])
                        ISF_matrix = calc_ISF_matrix(fai_matrix, sigema, miu, rou_up)
                        # print('ISF_matrix[rou][rou_up]',rou,rou_up,'sigema,miu=',sigema,miu)
                        eigenvalues, eigenvectors = np.linalg.eigh(ISF_matrix)
                        eigenvalues_int = [int(np.around(ee, 0)) for ee in eigenvalues]
                        # eigenvalues,eigenvectors = np.linalg.eig(ISF_matrix)
                        eigenvectors = eigenvectors.T  # 小贱人numpy非要按列写
                        # fai_matrix = cut_faimatrix_beta(fai_matrix) # hard code了为了补前面的锅，这是原本calc_fai_matrix得到的fai_matrix
                        beta_dict = {}
                        # same_lamada_dict, lamada_beta_counter_dict = calc_betas_max_dict(eigenvalues_int)
                        # # same_lamada_dict = {lamada: [number, [indexs]]}
                        # # lamada_beta_counter_dict = {lamada: 0}
                        # finish_schmitt_dict = {}
                        for l in range(len(eigenvalues)):
                            save_dict = {}
                            # eigenvalue_int = int(np.around(eigenvalues[l], 0))  # int化本征值，正确
                            eigenvalue_int = eigenvalues_int[l]
                            # if same_lamada_dict[eigenvalue_int][0] == 1:  # 一重根不用管beta,但要调整
                            #     eigenvector = eigenvectors[l]
                            # elif same_lamada_dict[eigenvalue_int][0] >= 2:  # 多重根先所有求和造第一个全对称，再挨个施密特正交化手续
                            #     eigenvector = schmitt_ly_and_renew_counter(eigenvalue_int, eigenvectors,
                            #                                                lamada_beta_counter_dict, same_lamada_dict,
                            #                                                finish_schmitt_dict)
                            # else:
                            #     print('sth about beta wrong!look the eigenvalues and the 3 dicts')
                            #     print(error)
                            # lamada_r_up = get_lamada_r(n_up, rou_up)
                            # lamada_r = eigenvalue_int + lamada_r_up
                            # # rou_new = from_lamada_get_r(lamada_r,n)
                            # # print('lamada_r,eigenvalues, rou_up, n',lamada_r,eigenvalues,rou_up, n)
                            # rou_new = from_lamada_rup_get_r(lamada_r, rou_up, n)
                            # beta = from_betadict_rounew_calc_beta(beta_dict, lamada_r)
                            # 绝对，相对相位判断，并再整形eigenvector
                            # TODO 按公式来说，此处整形的是eigenvector，而不是ISF_vector
                            if jugde_No_0_branchlaw(rou_new, rou_up) == 0:  # 绝对
                                if jugde_absolute_phase(eigenvector) == -1:  # 反了就改，对了不操作
                                    eigenvector = -1 * eigenvector
                                    # ISF_vector = np.sign(eigenvector) * eigenvector * eigenvector
                            else:  # 用相对相位 # 调用绝对相位一定要从存起来的ISF调用
                                rou_up_right = get_first_branchlaw(rou_new)
                                fai_right, ISF_right = get_right_ISF_first_no_0(sigema, miu, rou_new, beta,
                                                                                rou_up_right)
                                ISF_vector_right = np.sign(ISF_right) * np.sqrt(abs(ISF_right))
                                # print('fai_right=',fai_right)
                                m_now = get_m_now(rou_new, rou_up)
                                m_right, m_up_right = get_m_right(rou_new, rou_up, rou_up_right, m_up_now=1, m_up_up=1)
                                ISF_vector_now = calc_ISF_now(eigenvector, sigema, miu, rou_new, beta, fai_right,
                                                              rou_up_right, fai_matrix, rou_up, m_right, m_up_right,
                                                              m_now)
                                # TODO 更安全的if
                                #if np.around(abs(ISF_vector_right), 3) == np.around(abs(ISF_vector_now), 3):
                                if np.sign(ISF_vector_right) * np.sign(ISF_vector_now) == -1:
                                    eigenvector = -1 * eigenvector
                                # if np.sign(ISF_right)*np.sign(ISF_now) == -1:
                                # eigenvector = -1 * eigenvector
                                # ISF_vector = -1 * ISF_vector
                            ISF_vector = np.sign(eigenvector) * eigenvector * eigenvector
                            # save_dict['N'] = ISF_N_int
                            for k in range(len(fai_matrix)):
                                key = (tuple(fai_matrix[k][0]), tuple(fai_matrix[k][1]), fai_matrix[k][2])
                                # save_dict[key] = ISF_vector_int[k]
                                save_dict[key] = ISF_vector[k]
                            # print('save_ISF=', save_dict)
                            # print('eigenvalue=', eigenvalue_int, 'lamada_r_up', lamada_r_up, 'lamada_r', lamada_r)
                            # print('vector=', eigenvector)
                            # print('ISF_vector=', ISF_vector)
                            folder_path = _folder_path + 'rou' + str(rou_new) + '_beta' + str(beta) + '/'
                            file_path = _folder_path + 'rou' + str(rou_new) + '_beta' + str(beta) + '/' + "rou'" + str(
                                rou_up) + "_beta'" + str(1)
                            if os.path.exists(folder_path):  # 有就跳过
                                pass
                            else:
                                os.makedirs(folder_path)
                                print('=' * 10, '建立目录', folder_path, '=' * 10)
                            save_ISF_by_path(file_path, save_dict)
                            #
                            #
                            #
                            #
                            #
                            # TODO 不一致的参数有rou_new，beta，fai_matrix，save_dict
                            # 也就是说，也可以等ISF算完，一起算CGC
                            #
                            # 继承ISF结果，直接造CG系数
                            # ν:rou_new  ν':rou_up  [[σ',μ']]:fai_matrix  ISF:save_dict  β:beta  #β':beta_up
                            max_m_rou_up = get_Mm(Mm_path, rou_up, n_up)  # m'最大值
                            for m_rou_up in range(1, max_m_rou_up + 1):
                                m_rou = from_r_rup_mup_get_m(rou_new, rou_up, m_rou_up)
                                CG_file_name = _CG_folder_path + 'rou' + str(rou_new) + '_beta' + str(
                                    beta) + '_m' + str(
                                    m_rou) + '.pkl'
                                if os.path.exists(path):
                                    My_file = open(path, 'rb')
                                    CG_part = pickle.load(My_file)
                                    My_file.close()
                                else:
                                    CG_part = {'N': 1}
                                CG_preparation = copy.deepcopy(CG_part)
                                # νm:m_rou  ν'm':m_rou_up
                                for kk in range(len(fai_matrix)):
                                    fai_list = fai_matrix[kk]
                                    fai_tuple = (tuple(fai_matrix[kk][0]), tuple(fai_matrix[kk][1]), fai_matrix[kk][2])
                                    sigema_up_new = fai_list[0]
                                    miu_up_new = fai_list[1]
                                    # CG_up_new_path = CG_up_path+'sigema'+str(sigema_up_new)+'_'+'miu'+str(miu_up_new)+'/'+str(rou_up)+str(m_rou_up)+'.pkl'
                                    CG_up_new_path = CG_up_path + 'sigema' + str(sigema_up_new) + '_' + 'miu' + str(
                                        miu_up_new) + '/'  # 太长了，拼两次吧
                                    CG_up_new_path = CG_up_new_path + 'rou' + str(rou_up) + '_beta' + str(
                                        fai_matrix[kk][2]) + '_m' + str(m_rou_up) + '.pkl'
                                    CG_up_dict = get_CG(CG_up_new_path)
                                    CG_up_N = CG_up_dict['N']
                                    CG_up_dict.pop('N')
                                    # σ'_list:sigema_up_new  μ'_list:miu_up_new
                                    # σ'_tuple:fai_tuple[0]  μ'_tuple:fai_tuple[1]
                                    for m_sigema_up, m_miu_up in CG_up_dict:
                                        CG_up = CG_up_dict[(m_sigema_up, m_miu_up)]
                                        m_sigema = from_r_rup_mup_get_m(sigema, sigema_up_new, m_sigema_up)
                                        m_miu = from_r_rup_mup_get_m(miu, miu_up_new, m_miu_up)
                                        ISF = save_dict[(fai_tuple[0], fai_tuple[1], fai_tuple[2])]
                                        # TODO 先写上float版吧
                                        CG_float = ISF * CG_up / CG_up_N  # "用isf_square而不是isf倒是无所谓，因为它作用于整体"
                                        if (m_sigema, m_miu) in CG_preparation:
                                            CG_preparation[(m_sigema, m_miu)] = CG_preparation[(m_sigema, m_miu)] + CG_float  # "这里，用square乘完再加不太对吧？"
                                        else:
                                            CG_preparation[(m_sigema, m_miu)] = CG_float
                                # TODO 整形CG
                                CG_preparation.pop('N')
                                CG = CG_preparation
                                CG_N = 0
                                for ii in CG:
                                    CG_N += abs(CG[ii])
                                CG['N'] = CG_N
                                save_CG_by_path(CG_file_name, CG)
                                # print('save:CG_now=', CG, CG_N)

                    if not np.around(CG_N, 0) == 1:  # 检查最后CG的归一性  "这个地方要紧一点，毕竟0.51~1.49都能通过" "例如卡到np.around(CG_N, 1) 就收紧到0.950~1.049"
                        print('Wrong! and CG_N=', CG_N, 'CG=', CG)
                        print(error)
