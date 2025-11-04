import os
import numpy as np
import math
from algorithm import *


P_mass = 1.007277  # 质子相对质量
H2O_mass = 1.0078246*2 + 15.9949141  # 水分子相对质量

"""
计算precursor mh，不支持修饰
"""
def cal_precursor_mh(seq, aa2mass):
    mh = 0.0
    for aa in seq:
        mh += aa2mass[aa]
    mh += P_mass
    mh += H2O_mass
    return mh


"""
计算precursor mh，支持修饰
"""
def cal_precursor_mh_(modified_seq, aa2mass, unimod2mass_shift):
    """
    使用方法:
        element2mass = get_element_mass("infos/element.ini")
        aa2mass = get_aa_mass("infos/aa.ini", element2mass)
        precursor_mh = cal_precursor_mh("FYPDRPHQK", aa2mass)
    """
    mh = 0.0
    parts = modified_seq.split("(")
    for part in parts:
        items = part.split(")")
        if len(items) > 1:  # unimod, seq
            unimod, seq = items
            mh += unimod2mass_shift[unimod]
            for aa in seq:
                mh += aa2mass[aa]
        else:
            for aa in items[0]:
                mh += aa2mass[aa]

    mh += P_mass
    mh += H2O_mass
    return mh


"""
# 元素-分子量字典
get_element2mass("infos/element.ini")
"""
def get_element2mass(element_ini_path):
    f = open(element_ini_path, 'r')
    element2mass = dict()
    for line in f.readlines()[1:]:
        items = line.strip('|').split('|')
        name = items[0].split('=')[1]
        masses = items[1].strip(',').split(',')
        abundances = items[2].strip(',').split(',')
        element_mass = 0
        for mass, abundance in zip(masses, abundances):
            # average
            # element_mass += float(mass)*float(abundance)
            # mono
            element_mass += float(mass)
            break
        element2mass[name] = element_mass
    element2mass['P'] = 1.007277  # 质子P的质量
    f.close()
    return element2mass


"""
# 氨基酸残基-分子量字典
get_aa2mass("infos/aa.ini", element2mass)
"""
def get_aa2mass(aa_ini_path, element2mass):
    f = open(aa_ini_path, 'r')
    aa2mass = dict()
    for line in f.readlines()[2:]:  # 从第3行开始
        items = line.strip('|').split('|')
        name = items[0].split('=')[1]
        element_nums = items[1].strip(')').split(')')
        aa_mass = 0
        for element_num in element_nums:
            element_num_list = element_num.split('(')
            element = element_num_list[0]
            num = int(element_num_list[1])
            aa_mass += element2mass[element]*num
        aa2mass[name] = aa_mass
    f.close()
    return aa2mass

"""
# 氨基酸残基-分子量字典
get_mod2mass("infos/modification.ini")
"""
def get_mod2mass(mod_ini_path):
    f = open(mod_ini_path, "r")
    mod2mass = dict()
    next(f,None)
    for line in f:
        if line.startswith("name"):
            continue
        mod_name, info = line.split("=")
        mod2mass[mod_name] = float(info.split()[2])
    f.close()
    return mod2mass
        
    
def get_unimod2mass(mod_ini_path, element2mass):
    return {"UniMod:1": 42.010565,
            "UniMod:4": 57.021464,
            "UniMod:35": 15.994915}


"""
计算modified_seq与谱图的离子匹配情况: modification.ini
√ 考虑b,y离子
x 考虑失水, 失氨离子
√ 考虑[1,precursor_charge)的所有电荷数目
- peaks为[mz,intensity]对的列表，按mz从小到大排序
"""
def cal_by_matches(peaks, modified_seq, precursor_charge, aa2mass, mod2mass, topk=None, ppm=20):
    """
    使用方法:
        cal_by_matches(peaks, "PANITDSLSNR", "PANITDSLSNR", 3, topk=200)
        return: (matched, seq_ratio, b_match_cnt, b_ratio, y_match_cnt, y_ratio)
    """
    if topk is not None:
        intensities = [peak[1] for peak in peaks]
        intensity_th = quickselect(intensities, topk)
        if intensity_th != -1:
            _peaks = list()
            for mz, intensity in peaks:
                if intensity >= intensity_th:
                    _peaks.append([mz, intensity])
            peaks = _peaks
            
    def search_mz(peaks, target_mz, ppm=20):
        """
        如果找到匹配峰，则返回该峰的mz, intensity, error, 否则返回-1, -1, -1
        """
        left_mz, right_mz = (1 - ppm * 1e-6) * target_mz, (1 + ppm * 1e-6) * target_mz
        left, right = 0, len(peaks)  # 左闭右开区间
        while left < right:
            mid = (left + right) // 2
            if left_mz <= peaks[mid][0] <= right_mz:
                return peaks[mid][0], peaks[mid][1], (target_mz-peaks[mid][0])/peaks[mid][0]
            elif peaks[mid][0] < left_mz:
                left = mid + 1
            else:
                right = mid
        return -1,-1,-1

    seq, mods = modified_seq.split("\t")
    n = len(seq)

    # 计算每个位点处，氨基酸+修饰的质量
    r_mass = [0.0] * n
    for i in range(n):
        r_mass[i] = aa2mass[seq[i]]
    for mod in mods.split(";"):
        if mod=='':
            break
        mod_pos,mod_name = mod.split(",")
        mod_pos = int(mod_pos)
        if mod_pos==0:  # N-term
            r_mass[0] += mod2mass[mod_name]
        else:  # normal & C-term
            r_mass[mod_pos-1] += mod2mass[mod_name] 

    matched = list()  # 匹配情况
    b_mz = 0.0
    for i in range(n - 1):
        b_mz += r_mass[i]
        for charge in range(1,precursor_charge+1):
            matched_mz, matched_intensity, error = search_mz(peaks, b_mz/charge+P_mass, ppm=ppm)
            if matched_intensity >= 0:  # 有匹配峰
                matched.append((f"b{i+1}"+("+"*charge), matched_mz, matched_intensity, error))  # (ion_name,mz,intensity,error)
    y_mz = H2O_mass
    for i in range(n - 1, 0, -1):
        y_mz += r_mass[i]
        for charge in range(1,precursor_charge+1):
            matched_mz, matched_intensity, error = search_mz(peaks, y_mz/charge+P_mass, ppm=ppm)
            if matched_intensity >= 0:
                matched.append((f"y{n-i}"+("+"*charge), matched_mz, matched_intensity, error))

    return matched


"""
计算modified_seq与谱图的离子匹配情况: modification.ini
√ 考虑b,y离子
x 考虑失水, 失氨离子
√ 考虑[1,precursor_charge)的所有电荷数目
- peaks为[mz,intensity]对的列表，按mz从小到大排序
"""
def cal_by_matches_da(peaks, modified_seq, precursor_charge, aa2mass, mod2mass, topk=None, da=0.02):
    """
    使用方法:
        cal_by_matches(peaks, "PANITDSLSNR", "PANITDSLSNR", 3, topk=200)
        return: (matched, seq_ratio, b_match_cnt, b_ratio, y_match_cnt, y_ratio)
    """
    if topk is not None:
        intensities = [peak[1] for peak in peaks]
        intensity_th = quickselect(intensities, topk)
        if intensity_th != -1:
            _peaks = list()
            for mz, intensity in peaks:
                if intensity >= intensity_th:
                    _peaks.append([mz, intensity])
            peaks = _peaks
            
    def search_mz(peaks, target_mz, da=0.02):
        """
        如果找到匹配峰，则返回该峰的mz, intensity, error, 否则返回-1, -1, -1
        """
        left_mz, right_mz = target_mz-da, target_mz+da
        left, right = 0, len(peaks)  # 左闭右开区间
        while left < right:
            mid = (left + right) // 2
            if left_mz <= peaks[mid][0] <= right_mz:
                return peaks[mid][0], peaks[mid][1], (target_mz-peaks[mid][0])/peaks[mid][0]
            elif peaks[mid][0] < left_mz:
                left = mid + 1
            else:
                right = mid
        return -1,-1,-1

    seq, mods = modified_seq.split("\t")
    n = len(seq)

    # 计算每个位点处，氨基酸+修饰的质量
    r_mass = [0.0] * n
    for i in range(n):
        r_mass[i] = aa2mass[seq[i]]
    for mod in mods.split(";"):
        if mod=='':
            break
        mod_pos,mod_name = mod.split(",")
        mod_pos = int(mod_pos)
        if mod_pos==0:  # N-term
            r_mass[0] += mod2mass[mod_name]
        else:  # normal & C-term
            r_mass[mod_pos-1] += mod2mass[mod_name] 

    matched = list()  # 匹配情况
    b_mz = 0.0
    for i in range(n - 1):
        b_mz += r_mass[i]
        for charge in range(1,precursor_charge+1):
            matched_mz, matched_intensity, error = search_mz(peaks, b_mz/charge+P_mass, da=da)
            if matched_intensity >= 0:  # 有匹配峰
                matched.append((f"b{i+1}"+("+"*charge), matched_mz, matched_intensity, error))  # (ion_name,mz,intensity,error)
    y_mz = H2O_mass
    for i in range(n - 1, 0, -1):
        y_mz += r_mass[i]
        for charge in range(1,precursor_charge+1):
            matched_mz, matched_intensity, error = search_mz(peaks, y_mz/charge+P_mass, da=da)
            if matched_intensity >= 0:
                matched.append((f"y{n-i}"+("+"*charge), matched_mz, matched_intensity, error))

    return matched



"""
给定DIA-NN形式的修饰肽段序列，输出未修饰肽段序列
"""
def remove_unimod(modified_seq):
    seq = ""
    parts = modified_seq.split('(')
    for part in parts:
        seq += part.split(")")[-1]
    return seq


"""
获取pfind_idx -> peaks
"""
def get_idx2peaks(mgf_paths):
    idx2peaks = dict()  # pfind_idx->peaks
    for mgf_path in mgf_paths:
        print(mgf_path)
        fm = open(mgf_path, "r")
        while True:
            line = fm.readline()
            if not line:
                break
            if line.startswith("B"):  # BEGIN IONS
                pass
            elif line.startswith("C"):  # CHARGE
                charge = int(line.strip().split("=")[-1].split("+")[0])  # 3+
            elif line.startswith("T"):  # TITLE
                title = line.strip().split("=")[-1]
                file_id = shorten_name(title)
                idx = title.split(".")[-4]
                pfind_idx = file_id + "-" + idx
            elif line.startswith("P"):  # PEPMASS
                mz = float(line.strip().split("=")[-1])
            # elif line.startswith("E"):  # END IONS
            #     pass
            elif "A" <= line[0] <= "Z":
                continue
            else:  # PEAKS
                peaks = list()
                while True:
                    line = fm.readline()
                    if line.startswith("E"):
                        break
                    items = line.strip().split()
                    mz = float(items[0])
                    intensity = float(items[1])
                    peaks.append((mz, intensity))
                idx2peaks[pfind_idx] = peaks
        fm.close()
    return idx2peaks


"""
input: res_paths(raw_x.txt)
output: pfind_idx -> [line.split(),...]
"""
def get_pfind_raw_res(res_paths, idx=False):
    idx2pfind_raw_res = dict()
    for res_path in res_paths:
        print(res_path)
        f = open(res_path,"r")
        res = None
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith("b"):
                items = line.strip().split("\t")
                if idx:
                    es = items[1].split(".")
                    pfind_idx = es[1]+"-"+es[4]
                else:
                    pfind_idx = items[1]
                charge = int(items[2])
                mz = float(items[3])
                
                res = list()
            elif line.startswith("e"):
                idx2pfind_raw_res[pfind_idx] = res
            else:
                res.append(line.strip().split("\t"))
        f.close()
        
    return idx2pfind_raw_res


"""
input: res_paths(.qry.res or .sa)
output: pfind_idx -> [line.split(),...]
"""
def get_pfind_topk_res(res_paths, idx=False):
    idx2pfind_topk_res = dict()
    for res_path in res_paths:
        print(res_path)
        f = open(res_path,"r")
        res = None
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith("S"):
                if res is not None and len(res)>0:
                    idx2pfind_topk_res[pfind_idx] = res
                _, mz, charge = line.strip().split("\t")
                mz = float(mz)
                charge = int(charge)

                line = f.readline()
                if idx:
                    items = line.strip().split(".")
                    pfind_idx = items[1] + "-" + items[4]
                else:
                    pfind_idx = line.strip()

                res = list()
            else:
                res.append(line.strip().split("\t"))
        if res is not None and len(res)>0:
            idx2pfind_topk_res[pfind_idx] = res
        f.close()
        
    return idx2pfind_topk_res

"""
input: res_paths(.spectra)
"""
def get_pfind_res(res_paths, idx=False):
    idx2pfind_res = dict()
    for res_path in res_paths:
        print(res_path)
        f = open(res_path,"r")
        next(f, None)
        res = None
        while True:
            line = f.readline()
            if not line:
                break
            items = line.strip().split("\t")
            if idx:
                es = items[0].split(".")
                pfind_idx = es[1]+"-"+es[4]
            else:
                pfind_idx = items[0]
            idx2pfind_res[pfind_idx] = items[1:]
        f.close()
        
    return idx2pfind_res
    
    
# pfind的scan_no向DIA-NN映射
# shortened_name->ori_scan_nos->new_scan_nos
def get_scan_mapping(pfc_paths):
    pfind_ms2toms1 = dict()  # DIA-NN二级质谱scan号对应的一级质谱scan号: pFind_ms2 -> pFind_ms1
    ms2_diann2pfind = dict()  # DIA-NN到pFind二级质谱scan号的映射: DIA-NN_ms2 -> pFind_ms2
    # diann_ms2toms1 = dict()  # DIA-NN二级质谱scan号对应的一级质谱scan号: DIA-NN_ms2 -> pFind_ms1
    # ms2_pfind2diann = dict()  # pFind到DIA-NN二级质谱scan号的映射: pFind_ms2 -> DIA-NN_ms2
    # pfind_ms1s = dict()  # shortened_name -> ms1s: top50_01 -> [ms1,...]
    for pfc_path in pfc_paths:
        shortened_name = shorten_name(pfc_path)
        print(shortened_name)
        
        # p2d, dms2toms1 = dict(), dict()
        d2p, pms2toms1 = dict(), dict()
        # ms1s = list()
        
        f = open(pfc_path, 'r')
        diann_ms2_no = 0  # DIA-NN的MS2从0开始编号
        for line in f.readlines()[1:]:
            items = line.strip().split('\t')
            if items[3] == "MS1":
                # ms1s.append(items[0])
                continue
            pfind_ms2_no = items[0]
            pfind_ms1_no = items[1]
            d2p[str(diann_ms2_no)] = pfind_ms2_no
            pms2toms1[items[0]] = pfind_ms1_no
            # p2d[pfind_ms2_no] = str(diann_ms2_no)
            # dms2toms1[str(diann_ms2_no)] = pfind_ms1_no
            diann_ms2_no += 1
        ms2_diann2pfind[shortened_name] = d2p
        pfind_ms2toms1[shortened_name] = pms2toms1
        # ms2_pfind2diann[shortened_name] = p2d
        # diann_ms2toms1[shortened_name] = dms2toms1
        # raw2ms1s[shortened_name] = ms1s
        f.close()
        
    return ms2_diann2pfind, pfind_ms2toms1


# 缩短raw文件名
def shorten_name(raw_name):
    """
    :param raw_name: 20230108_AST_Neo1_DDA_UHG_HeLa_200ng_2th2p5ms_top50_01_20230808155512.raw
    :return: top50_01
    """
    items = raw_name.strip('.raw').split('ms_')[-1]
    items = items.split('_')
    return items[0]+'_'+items[1].split('.')[0]


# 将modification.ini中的修饰映射至Unimod_id
def get_mod2unimod_id(mod_ini_path, unimod_path):
    funi = open(unimod_path, "r")
    name2id = dict()
    for line in funi.readlines():
        if line.startswith("id"):
            unimod_id = line.strip().split(": ")[1]
        if line.startswith("name"):
            unimod_name = line.strip().split(": ")[1].lower().replace(":", "_")  # 注意转小写，替换:为_
            name2id[unimod_name] = unimod_id

    fini = open(mod_ini_path, "r")
    mod2unimod_id = dict()
    for line in fini.readlines():
        if line.startswith("name"):
            mod_name = line.strip().split()[0].split('=')[1]
            unimod_name = "".join(mod_name.split('[')[:-1]).lower()                    # 麻烦死了T_T
            if unimod_name in name2id:  # 注意，unimod中可能没有modification.ini的某些修饰
                mod2unimod_id[mod_name] = name2id[unimod_name].replace('UNIMOD','UniMod')  # 神经，这都不一样
            else:
                mod2unimod_id[mod_name] = "UniMod:-1"  # 此时标记为unknown就行，反正DIA-NN没有几个修饰
    return mod2unimod_id


# 将pfind的肽段+修饰结果转换为DIA-NN的格式
def mod_pfind2diann(seq, mod, mod2unimod_id):
    """
    :param seq: pfind预测的肽段
    :param mod: pfind预测的修饰
    :param mod2unimod: modification.ini中修饰名称与UniMod对应id(accession)的映射
    """
    if mod == "":
        return seq

    pfind_modified_seq = ""
    last_pos = 0
    for pos_name in mod.strip(';').split(';'):
        pos, mod_name = pos_name.split(',')
        pos = int(pos)
        unimod_id = mod2unimod_id[mod_name]
        for i in range(last_pos, min(pos,len(seq))):  # 注意可能有C端修饰，不过DIA-NN中没有此修饰，所以我暂时把修饰直接加在末尾
            pfind_modified_seq += seq[i]
        pfind_modified_seq += ('('+unimod_id+')')
        last_pos = pos
    pfind_modified_seq += seq[last_pos:]
    return pfind_modified_seq


def get_pfind_result(pfind_res_path, mod2unimod_id, pep_level_fdr=0.01, scan2diann=False, scan_pfind2diann=None):
    fr = open(pfind_res_path, "r")
    scan2info = dict()
    for line in fr.readlines()[1:]:
        items = line.strip().split("\t")
        q_value = float(items[4])
        if len(items) <= 5 or q_value > pep_level_fdr:
            break
        raw_name = items[0].split('\\')[-1]
        shortened_name = shorten_name(raw_name)
        scan_no = items[1] if not scan2diann else scan_pfind2diann[shortened_name][items[1]]
        idx = shortened_name + '-' + scan_no

        exp_mh = float(items[2])
        charge = int(items[3])
        seq = items[5]
        cal_mh = float(items[6])
        mod = items[10]
        modified_seq = mod_pfind2diann(seq, mod, mod2unimod_id)
        protein_names = items[12].strip('/').split('/')
        # seqs.add(seq+mod)
        is_decoy = (0 if items[15]=='target' else 1)

        # 注意：很有可能一个谱图中包含多个结果(混合谱)
        if idx not in scan2info:
            scan2info[idx] = [[seq,modified_seq,charge,is_decoy,cal_mh,exp_mh,q_value,protein_names]]
        else:
            # print('mixture spectrum:', idx)
            scan2info[idx].append([seq,modified_seq,charge,is_decoy,cal_mh,exp_mh,q_value,protein_names])
    # print(len(seqs))
    return scan2info


def get_diann_result(diann_report_path, diann_reportlib_path, scan2pfind=True, scan_diann2pfind=None):
    seq2info = dict()  # pep_seq->[肽段理论质荷比(cal_mz),is_decoy]
    with open(diann_reportlib_path, "r") as frl:
        next(frl, None)  # 跳过第一行
        for line in frl:
            items = line.strip().split("\t")
            pep_seq = items[9]
            modified_seq = items[17]
            if modified_seq not in seq2info:
                cal_mz = float(items[1])
                is_decoy = int(items[8])
                charge = int(items[19])
                cal_mh = charge*cal_mz - (charge-1)*1.007277
                # seq2info[pep_seq] = [cal_mz,is_decoy]
                seq2info[modified_seq] = [cal_mh,is_decoy]

    # seqs = set()
    scan2info = dict()  # raw_name+scan_no->[modified_seq,cal_mz,charge,is_decoy,q_value]
    with open(diann_report_path, "r") as fr:
        next(fr, None)  # 跳过第一行
        for line in fr:
            items = line.strip().split("\t")
            protein_ids_str = items[3]  # protein.ids
            seq = items[14]

            raw_name = items[0].split('\\')[-1]
            shortened_name = shorten_name(raw_name)
            scan_no = items[53] if not scan2pfind else scan_diann2pfind[shortened_name][items[53]]
            idx = shortened_name + '-' + scan_no

            modified_seq = items[13]
            # seqs.add(modified_seq)
            charge = int(items[16])
            q_value = float(items[17])

            # cal_mz = seq2info[seq][0]
            # is_decoy = seq2info[seq][1]
            cal_mh = seq2info[modified_seq][0]
            cal_mz = (cal_mh+(charge-1)*1.007277)/charge
            is_decoy = seq2info[modified_seq][1]

            # 注意：很有可能一个谱图中包含多个结果(混合谱)
            if idx not in scan2info:
                scan2info[idx] = [[seq,modified_seq,charge,is_decoy,cal_mz,q_value,protein_ids_str]]
            else:
                scan2info[idx].append([seq,modified_seq,charge,is_decoy,cal_mz,q_value,protein_ids_str])

        return scan2info


# 获取蛋白质名称->seq的映射
def get_pro2seq(fasta_path):
    pro2seq = dict()
    seq, pro_name = "", None
    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith('>'):
                if pro_name:
                    pro2seq[pro_name] = seq
                pro_name = line[1:].strip().split()[0]  # 简单地取第1个元素为protein name
                seq = ""
            else:
                seq += line.strip()
    if seq != "":
        pro2seq[pro_name] = seq
    return pro2seq


def get_uniprot_id2pro_name(fasta_path):
    uniprot_id2pro_name = dict()
    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith('>'):
                pro_name = line[1:].strip().split()[0]  # 简单地取第1个元素为protein name
                if pro_name.startswith("sp"):
                    uniprot_id = pro_name.split('|')[1]
                elif pro_name.startswith("CON"):
                    uniprot_id = pro_name
                else:
                    print("error!")
                uniprot_id2pro_name[uniprot_id] = pro_name
    return uniprot_id2pro_name


# 根据序列seq，找到其对应的蛋白质: replace=True表示对seq的I/L进行替换
def get_seq2pro(seq, pro2seq, replace=True):
    protein_names = list()
    seqs = [seq]
    if replace:
        for i in range(len(seqs)):
            if seqs[i] == "I":
                seq_temp = seq.copy()
                seq_temp[i] = "L"
                seqs.append(seq_temp)
            elif seqs[i] == "L":
                seq_temp = seq.copy()
                seq_temp[i] = "I"
                seqs.append(seq_temp)

    for seq in seqs:
        for pro,pro_seq in pro2seq.items():
            if seq in pro_seq:
                protein_names.append(pro)
    return protein_names


# # 将protien_ids_str解析为protein_names(list), 屎一样的DIA-NN竟然拿';'作蛋白质的分隔符(很多蛋白质名称也包含';'啊!!!)
# def parse(protein_ids_str):
#     protein_ids = protein_ids_str.split(';')
#     protein_names = list()
#     pre = ""
#     for protein_id in protein_ids:
#         protein_id = pre + protein_id
#         if protein_id in uniprot_id2pro_name:
#             protein_names.append(protein_id)
#             pre = ""
#         else:
#             pre = protein_id + ";"
#     if pre != "" or len(protein_names) == 0:
#         print("protein_ids_str parse error!", protein_ids_str)
#     return protein_names


def parse(protein_ids_str, uniprot_id2pro_name):
    protein_names = list()
    semis_pos = [-1]
    for pos, s in enumerate(protein_ids_str):
        if s == ";":
            semis_pos.append(pos)
    semis_pos.append(len(protein_ids_str))  # [0,...,n]

    p1, p2 = 0, len(semis_pos) - 1
    while p1 < p2:
        pro_id = protein_ids_str[semis_pos[p1] + 1:semis_pos[p2]]
        if pro_id in uniprot_id2pro_name:
            protein_names.append(uniprot_id2pro_name[pro_id])
            p1 = p2
            p2 = len(semis_pos) - 1
        else:
            p2 -= 1

    if p1 == p2 == len(semis_pos) - 1 and len(protein_names) != 0:
        return protein_names
    else:
        print("protein_ids_str parse error!", protein_ids_str)


"""
解析pParse导出的母离子: 
scan2precursors_info: pfind_idx -> [(mz,charge)...]
"""
def get_scan2precursors_info(csv_paths):
    scan2precursors_info = dict()
    for csv_path in csv_paths:
        shortened_name = shorten_name(csv_path)
        with open(csv_path, "r") as f:
            next(f, None)
            for line in f:
                items = line.strip(", ").split(", ")

                scan_no = items[2]

                # diann_scan_no = scan_pfind2diann[shortened_name][scan_no]
                pfind_idx = shortened_name + '-' + scan_no

                scan2precursors_info[pfind_idx] = list()
                precursor_num = int(items[5])
                for i in range(precursor_num):  # 每三个一组，分别为idx,mz,charge
                    mz, charge = float(items[5 + i * 3 + 2]), int(items[5 + i * 3 + 3])
                    scan2precursors_info[pfind_idx].append((mz, charge))
    return scan2precursors_info


def get_pro2category(category_dir):
    """
    获取蛋白质类别
    category_dir: 该文件夹下只允许存放类别文件(以".txt"结尾)，每个文件中存放相同类型的蛋白质名称，
                  注意蛋白质名称要与.pac文件中的一致
    """
    pro2category = dict()
    for file_name in os.listdir(category_dir):
        category = file_name.strip(".txt")  # 蛋白质类型
        file_path = os.path.join(category_dir, file_name)
        with open(file_path, "r") as f:
            for line in f:
                protein_name = line.strip()
                pro2category[protein_name] = category
                # decoy蛋白质
                rev_protein_name = "REV_" + protein_name
                pro2category[rev_protein_name] = "decoy"
    return pro2category


def report_pfind(res_path, pro2category):
    proteins = set()
    category2cnt = dict()
    for key in set(pro2category.values()):
        category2cnt[key] = 0

    f = open(res_path, "r")
    for line in f.readlines()[1:]:
        items = line.strip().split('\t')
        if not line.startswith('\t'):
            # flag = False
            if len(items) > 1:
                pro_name = items[1]
            else:
                break  # fdr=0.01
        elif line.startswith("	SubSet"):
            continue
        elif line.startswith("	SameSet"):
            # if items[1].endswith("_HUMAN") and pro_name.endswith("_MOUSE"):
            #     if not flag:
            #         category2cnt[pro2category[pro_name]] -= 1
            #         flag = True
            #     category2cnt[pro2category[items[1]]] += 1
            continue
        else:
            continue

        if pro_name not in proteins:
            proteins.add(pro_name)
            category2cnt[pro2category[pro_name]] += 1

    for category, cnt in category2cnt.items():
        print(category, cnt)
    return proteins


def report_my_res(res_path, pro2category):
    proteins = set()
    category2cnt = dict()
    for key in set(pro2category.values()):
        category2cnt[key] = 0

    f = open(res_path, "r")
    for line in f.readlines():
        pro_name = line.strip()

        if not pro_name in proteins:
            proteins.add(pro_name)
            category2cnt[pro2category[pro_name]] += 1

    for category, cnt in category2cnt.items():
        print(category, cnt)
    return proteins