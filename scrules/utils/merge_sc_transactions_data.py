# -*-coding:utf-8-*-
"""
@Author: Phantom
@Time:2023/9/5 17:48
@Email: 2909981736@qq.com
"""
import csv


def merge_files(sc_transactions_data_files, savepath):
    with open(savepath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for file in sc_transactions_data_files:
            with open(file, 'r') as infile:
                reader = csv.reader(infile)
                for row in reader:
                    writer.writerow(row)
    print("Merged successfully")
