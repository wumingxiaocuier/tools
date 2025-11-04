import os
import random


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


def quickselect(arr, k):
    """
    Find the k-th largest element in an unsorted list.
    :param arr: List of integers
    :param k: The k-th largest element to find
    :return: The k-th largest element
    """
    if not 1 <= k <= len(arr):
        return -1

    def partition(left, right, pivot_index):
        pivot = arr[pivot_index]
        arr[pivot_index], arr[right] = arr[right], arr[pivot_index]
        store_index = left
        for i in range(left, right):
            if arr[i] > pivot:
                arr[i], arr[store_index] = arr[store_index], arr[i]
                store_index += 1
        arr[right], arr[store_index] = arr[store_index], arr[right]
        return store_index

    def select(left, right, k_smallest):
        if left == right:
            return arr[left]
        pivot_index = random.randint(left, right)
        pivot_index = partition(left, right, pivot_index)
        if k_smallest == pivot_index:
            return arr[k_smallest]
        elif k_smallest < pivot_index:
            return select(left, pivot_index - 1, k_smallest)
        else:
            return select(pivot_index + 1, right, k_smallest)

    return select(0, len(arr) - 1, k - 1)




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