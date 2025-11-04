import os


def get_pi(p):
    pi = [-1] * len(p)
    p1, p2 = 0, -1  # p1用来遍历字符串，p2用作匹配
    while p1 < len(p)-1:
        if p2 == -1 or p[p1] == p[p2]:
            p1 += 1
            p2 += 1
            pi[p1] = p2
        else:
            p2 = pi[p2]
    return pi


def KMP(s, p):
    """
    查找p在s中出现的第一个位置下标
    """
    ps, pp = 0, 0
    lens, lenp = len(s), len(p)
    pi = get_pi(p)
    while ps < lens and pp < lenp:
        if pp == -1 or s[ps] == p[pp]:
            ps += 1
            pp += 1
        else:
            pp = pi[pp]
    if pp == lenp:
        return ps - pp
    else:
        return -1


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


# 大顶堆
class MaxHeap(object):
    def __init__(self, items=None):
        if items is None:
            self.items = list()
            self.size = 0
        else:
            self.items = items
            self.size = len(self.items)
        self._build()  # 建堆

    def _adjust_down(self, parent):
        son = parent*2 + 1
        temp = self.items[parent]
        while son < self.size:
            if son+1 < self.size and self.items[son+1] > self.items[son]:
                son += 1
            if temp < self.items[son]:
                self.items[parent] = self.items[son]
                parent = son
                son = parent*2 + 1
            else:
                break
        self.items[parent] = temp

    def _adjust_up(self, son):
        parent = (son-1) // 2
        temp = self.items[son]
        while parent >= 0:
            if self.items[parent] < temp:
                self.items[son] = self.items[parent]
                son = parent
                parent = (son-1) // 2
            else:
                break
        self.items[son] = temp

    def _build(self):
        """
        建堆
        """
        for i in reversed(range(self.size//2)):
            self._adjust_down(i)

    def push(self, item):
        self.items.append(item)
        self.size += 1
        self._adjust_up(self.size - 1)

    def pop(self):
        assert len(self.items) > 0, "Cannot pop item! The MaxHeap is empty!"
        self.items[0], self.items[self.size-1] = self.items[self.size-1], self.items[0]
        self.size -= 1
        self._adjust_down(0)

        return self.items.pop()


class MaxHeapForProteinInference(MaxHeap):
    """
    继承自MaxHeap，增加修改功能
    - items为Protein对象列表，每当某个Protein对象修改后，需要重新计算其在堆中的位置(需要保存一个protein_name->id的映射)
    """
    def __init__(self, items=None):
        self.protein_name2id = dict()
        super().__init__(items)
        for idx, protein in enumerate(self.items):
            self.protein_name2id[protein.name] = idx

    # 重写_adjust_down()
    def _adjust_down(self, parent):
        son = parent*2 + 1
        temp = self.items[parent]
        while son < self.size:
            if son+1 < self.size and self.items[son+1] > self.items[son]:
                son += 1
            if temp < self.items[son]:
                self.items[parent] = self.items[son]
                self.protein_name2id[self.items[son].name] = parent
                parent = son
                son = parent*2 + 1
            else:
                break
        self.items[parent] = temp
        self.protein_name2id[temp.name] = parent

    # 重写_adjust_up()
    def _adjust_up(self, son):
        parent = (son-1) // 2
        temp = self.items[son]
        while parent >= 0:
            if self.items[parent] < temp:
                self.items[son] = self.items[parent]
                self.protein_name2id[self.items[parent].name] = son
                son = parent
                parent = (son-1) // 2
            else:
                break
        self.items[son] = temp
        self.protein_name2id[temp.name] = son

    def adjust_protein(self, protein_name):
        """
        对于修改过的蛋白质，重新寻找其在堆中的位置
        注意，这个修改指的是"减小"，因此才选择__adjust_down(); 若是增大，则需要__adjust_up()
        """
        self._adjust_down(self.protein_name2id[protein_name])