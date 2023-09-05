# -*-coding:utf-8-*-
"""
@Author: Phantom
@Time:2023/9/5 17:50
@Email: 2909981736@qq.com
"""
import pandas as pd
from sqlalchemy import create_engine


def csv2sql(rules_file, username, password, ip, database, table_name):

    df = pd.read_csv(rules_file)
    if not df.isna().any().any():
        df.replace([float('inf'), -float('inf')], [9999, -9999], inplace=True)
        engine = create_engine(f'mysql://{username}:{password}@{ip}/{database}')
        df.to_sql(table_name, engine, if_exists='replace', index=False)
        print(table_name, 'ok!')

    else:
        print('Error: There are NaN values in the data frame and cannot be saved as sql')

