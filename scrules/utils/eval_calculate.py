# -*-coding:utf-8-*-
"""
@Author: Phantom
@Time:2023/9/4 17:51
@Email: 2909981736@qq.com
"""
import numpy as np


def calculate_conviction(support_B, confidence_AB):
    if confidence_AB == 1.0:
        conviction = float('inf')
    else:
        conviction = (1 - support_B) / (1 - confidence_AB)
    return conviction


def calculate_leverage(support_A, support_B, confidence_AB):
    support_AB = confidence_AB * support_A
    leverage = support_AB - (support_A * support_B)
    return leverage


def calculate_zhangs_metric(support_A, support_B, confidence_AB):
    support_AB = confidence_AB * support_A
    numerator = support_AB - (support_A * support_B)
    denominator = max(support_AB * (1 - support_A), support_A * (support_B - support_AB))
    if denominator == 0.0:
        if numerator == 0.0:
            zhangs_metric = 0.0
        else:
            zhangs_metric = numerator / abs(numerator) * float('inf')
    else:
        zhangs_metric = numerator / denominator
    return zhangs_metric


def calculate_chi_square(support_A, support_B, confidence_AB, total_transactions):
    support_AB = confidence_AB * support_A
    expected_AB = support_A * support_B * total_transactions
    observed_AB = support_AB * total_transactions

    chi_square = (observed_AB - expected_AB) ** 2 / expected_AB
    return chi_square


def calculate_eval_vectorization(rule_df):
    support_A = rule_df['antecedents_support']
    support_B = rule_df['consequents_support']
    confidence_AB = rule_df['confidence']
    support_AB = confidence_AB * support_A
    total_transactions = rule_df['total_transactions']
    expected_AB = support_A * support_B * total_transactions
    observed_AB = support_AB * total_transactions

    rule_df['conviction'] = (1 - support_B) / (1 - confidence_AB)
    rule_df.loc[rule_df['confidence'] == 1.0, 'conviction'] = float('inf')

    rule_df['leverage'] = support_AB - (support_A * support_B)

    rule_df['zhangs_metric'] = (support_AB - (support_A * support_B)) / (
        np.maximum(support_AB * (1 - support_A), support_A * (support_B - support_AB)))

    rule_df.loc[(rule_df['confidence'] * rule_df['antecedents_support'] - (rule_df['antecedents_support'] *
                                                                           rule_df['consequents_support']) > 0) & (
                        np.maximum(rule_df['confidence'] * rule_df['antecedents_support'] *
                                   (1 - rule_df['antecedents_support']),
                                   rule_df['antecedents_support'] * (rule_df['consequents_support'] -
                                                                     rule_df['confidence'] * rule_df[
                                                                         'antecedents_support'])) == 0), 'zhangs_metric'] = float(
        'inf')

    rule_df.loc[(rule_df['confidence'] * rule_df['antecedents_support'] - (rule_df['antecedents_support'] *
                                                                           rule_df['consequents_support']) < 0) & (
                        np.maximum(rule_df['confidence'] * rule_df['antecedents_support'] *
                                   (1 - rule_df['antecedents_support']),
                                   rule_df['antecedents_support'] * (rule_df['consequents_support'] -
                                                                     rule_df['confidence'] * rule_df[
                                                                         'antecedents_support'])) == 0), 'zhangs_metric'] = -float(
        'inf')

    rule_df.loc[(rule_df['confidence'] * rule_df['antecedents_support'] - (rule_df['antecedents_support'] *
                                                                           rule_df['consequents_support']) == 0) & (
                        np.maximum(rule_df['confidence'] * rule_df['antecedents_support'] *
                                   (1 - rule_df['antecedents_support']),
                                   rule_df['antecedents_support'] * (rule_df['consequents_support'] -
                                                                     rule_df['confidence'] * rule_df[
                                                                         'antecedents_support'])) == 0), 'zhangs_metric'] = 0.0

    rule_df['chi_square'] = (observed_AB - expected_AB) ** 2 / expected_AB

    return rule_df
