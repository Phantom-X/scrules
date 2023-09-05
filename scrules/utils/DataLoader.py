# -*-coding:utf-8-*-
"""
@Author: Phantom
@Time:2023/9/4 16:30
@Email: 2909981736@qq.com
"""
import csv


def load_sc_transactions_data(path):
    """
    Import single-cell data from a file that looks like this:
    item1, item2, item3,...
    item2,item6,item4,item5,...
    item3, item4, item5,...
    ...

    Note:
    - The file should be in comma-separated values (CSV) format.

    :param path:file path
    :return:sc_transactions_data:A two-dimensional list

    Example:
    >>> load_sc_transactions_data('sc_transactions.csv')
    [['item1', 'item2', 'item3'], ['item2', 'item6', 'item4','item5'], ...]
    """
    sc_transactions_data = []
    with open(path, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            row = list(set(row))
            row.sort()
            sc_transactions_data.append(row)
    return sc_transactions_data
