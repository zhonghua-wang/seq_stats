from collections import Counter
from typing import List, Tuple, Dict
import re

import pandas as pd
from pydantic import BaseModel
from tqdm import tqdm


# check result model
class CheckResult(BaseModel):
    start: int
    end: int
    type: str
    problem: str

    def add_base(self, base: int):
        return CheckResult(
            start=self.start + base, end=self.end + base,
            type=self.type, problem=self.problem
        )


# calculate the GC content of the specified sequence
def gc_content(seq: str) -> float:
    return (seq.count('G') + seq.count('C')) / len(seq)


# verify the sequence is valid, which contains only 'A', 'T', 'C', 'G'
def is_valid(seq: str) -> bool:
    return re.match(r'^[ATCG]+$', seq, re.IGNORECASE) is not None


# find the regions that GC contents < 40% or > 60% in the specified sequence
def gc_check(seq: str, window_size: int = 15, gc_threshold=(0.4, 0.6)) -> List[CheckResult]:
    lower, upper = gc_threshold
    gc_regions = []
    if len(seq) < window_size:
        # print('Invalid sequence')
        return []
    end = window_size
    for j in range(window_size, len(seq) + 1):
        # print(seq[start:j])
        if gc_content(seq[0:j]) < lower or gc_content(seq[0:j]) > upper:
            end = j
    gc_content_ = gc_content(seq[0:end])
    if gc_content_ < lower:
        gc_regions.append(CheckResult(start=0, end=end, type='GC含量', problem=f'GC含量 < {lower}'))
        gc_regions.extend([x.add_base(end) for x in gc_check(seq[end:], window_size, gc_threshold)])
    elif gc_content_ > upper:
        gc_regions.append(CheckResult(start=0, end=end, type='GC含量', problem=f'GC含量 > {upper}'))
        gc_regions.extend([x.add_base(end) for x in gc_check(seq[end:], window_size, gc_threshold)])
    else:
        gc_regions.extend([x.add_base(1) for x in gc_check(seq[1:], window_size, gc_threshold)])
    return gc_regions


# find the regions that contain consecutive nucleotides in the specified sequence
def consecutive_check(seq: str, window_size: int = 4, noise: float = 0.0) -> List[CheckResult]:
    if len(seq) < window_size:
        # print('Invalid sequence')
        return []
    consecutive_regions = []
    end = window_size
    for j in range(window_size, len(seq) + 1):
        # print(seq[start:j])
        mc = Counter(seq[0:j]).most_common()[0]
        if mc[1] / j >= 1 - noise:
            end = j
    mc = Counter(seq[0:end]).most_common()[0]
    if mc[1] / end >= 1 - noise:
        consecutive_regions.append(
            CheckResult(start=0, end=end, type='连续重复碱基', problem=f'连续碱基: {mc[0]}')
        )
        consecutive_regions.extend([x.add_base(end) for x in consecutive_check(seq[end:], window_size, noise)])
    else:
        consecutive_regions.extend([x.add_base(1) for x in consecutive_check(seq[1:], window_size, noise)])
    return consecutive_regions


# consecutive repeat unit check
# def repeat_unit_check_rex(seq: str, repeat_unit: str, repeat_count: int) -> List[CheckResult]:
#     # use regex to find the repeat unit
#     pattern = f"(AGC){{{repeat_count},}}"  # 匹配 AGC 重复 3 次以上的模式
#     match_list = []
#     for match in re.finditer(pattern, seq):
#         match_list.append(
#             CheckResult(start=match.start(), end=match.end(), type='unit_repeat', problem=f'repeat unit: {repeat_unit}')
#         )
#     return match_list


# consecutive repeat unit check with noise tolerance
def repeat_unit_check(seq: str, repeat_unit: str, repeat_count: int, noise: float = 0.0) -> List[CheckResult]:
    repeat_unit_list = [x.start() for x in re.finditer(repeat_unit, seq)]
    return _repeat_unit_check(seq, repeat_unit, repeat_unit_list, repeat_count, noise)


def _repeat_unit_check(
        seq: str, repeat_unit: str,
        repeat_unit_list: List[int],
        repeat_count: int, noise: float
) -> List[CheckResult]:
    if len(repeat_unit_list) < repeat_count:
        return []
    repeat_regions = []
    end = repeat_count - 1
    for j in range(repeat_count, len(repeat_unit_list)):
        no_repeat_nt_counts = repeat_unit_list[j] - repeat_unit_list[0] - j * len(repeat_unit)
        if no_repeat_nt_counts / (j + 1) <= noise:
            end = j
    no_repeat_nt_counts = repeat_unit_list[end] - repeat_unit_list[0] - end * len(repeat_unit)
    if no_repeat_nt_counts / (end + 1) <= noise:
        repeat_regions.append(
            CheckResult(
                start=repeat_unit_list[0], end=repeat_unit_list[end] + len(repeat_unit),
                type='连续重复序列', problem=f'重复单元: {repeat_unit}'
            )
        )
        repeat_regions.extend(
            [
                x
                for x in _repeat_unit_check(seq, repeat_unit, repeat_unit_list[end + 1:], repeat_count, noise)
            ]
        )
    else:
        repeat_regions.extend(
            [
                x
                for x in _repeat_unit_check(seq, repeat_unit, repeat_unit_list[1:], repeat_count, noise)
            ]
        )
    return repeat_regions


def _seq_cmp(seq1: str, seq2: str, noise: float):
    if noise == 0:
        return seq1 == seq2
    eq_pos = sum(map(lambda x: x[0] == x[1], (zip(seq1, seq2))))
    return eq_pos / len(seq1) >= 1 - noise


def long_distance_fixed_segment_repeats_check(dna_sequence: str, repeat_length: int, noise: float = 0) -> List[
    CheckResult]:
    """
    在DNA序列中查找多次出现的相同片段序列组以及它们的索引位置。

    参数：
    dna_sequence (str): DNA序列字符串。
    repeat_length (int): 相同片段序列的长度。

    返回：
    CheckResult列表，每个CheckResult对象包含一个相同片段序列组的信息。
    """
    repeats = {}
    for i in range(len(dna_sequence) - repeat_length + 1):
        repeat = dna_sequence[i:i + repeat_length]
        flag = False
        for k, v in repeats.items():
            # 如果重复片段与已有的重复片段相似，则将其添加到已有的重复片段组中
            if _seq_cmp(repeat, k, noise):
                flag = True
                pos_list = repeats[k]
                # append to repeats if the repeat is not adjacent to the previous one
                if i - pos_list[-1] - repeat_length > 1:
                    repeats[k].append(i)
                break
        if not flag:
            repeats[repeat] = [i]

    # 只保留重复出现的片段序列
    result_list = []
    for repeat, indexes in repeats.items():
        if len(indexes) > 1:
            for index in indexes:
                result_list.append(
                    CheckResult(
                        start=index, end=index + repeat_length,
                        type='非连续重复序列', problem=f'重复单元: {repeat}'
                    )
                )
    return result_list


def long_distance_segment_repeats_check(seq: str, min_length: int = 5, noise: float = 0) -> List[CheckResult]:
    result_list = []
    if len(seq) < min_length * 2:
        return result_list
    for length in range(int(len(seq) / 2), min_length - 1, -1):
        for result in long_distance_fixed_segment_repeats_check(seq, length, noise):
            if any(
                    [result.start <= x.start <= result.end or result.start <= x.end <= result.end for x in result_list]
            ):
                continue
            result_list.append(result)
    return result_list


def check_seq(
        seq: str,
        gc_min: float, gc_max: float, gc_window_size: int,
        consecutive_window_size: int, consecutive_noise: float,
        repeat_unit: str, repeat_count: int, repeat_noise: float,
        long_distance_repeats_min_length: int, long_distance_repeats_noise: float
) -> List[CheckResult]:
    # uppercase the sequence
    seq = seq.upper()
    # check if the sequence is valid
    if not is_valid(seq):
        raise ValueError('Invalid sequence.')
    result_list = []
    result_list.extend(gc_check(seq, gc_window_size, (gc_min, gc_max)))
    result_list.extend(consecutive_check(seq, consecutive_window_size, consecutive_noise))
    result_list.extend(repeat_unit_check(seq, repeat_unit, repeat_count, repeat_noise))
    result_list.extend(
        long_distance_segment_repeats_check(seq, long_distance_repeats_min_length, long_distance_repeats_noise)
    )
    return result_list


def check_fasta(
        fasta_str: str,
        gc_min: float, gc_max: float, gc_window_size: int,
        consecutive_window_size: int, consecutive_noise: float,
        repeat_unit: str, repeat_count: int, repeat_noise: float,
        long_distance_repeats_min_length: int, long_distance_repeats_noise: float
) -> Dict[str, List[CheckResult]]:
    seq_name = ''
    result_dict = {}
    for line in fasta_str.split('\n'):
        if line.startswith('>'):
            seq_name = line[1:].strip()
            # print(seq_name)
            continue
        else:
            seq = line.strip()
            if seq == '':
                continue
            if not is_valid(seq):
                print('Error (Invalid sequence): ', seq_name, seq)
            result_dict[seq_name] = []
            print(f'正在分析序列{seq_name} ...')
            # 把序列分成320bp一段，分别检查
            for i in range(0, len(seq), 320):
                # print(i)
                _results = check_seq(
                        seq[i:i + 320], gc_min, gc_max, gc_window_size, consecutive_window_size, consecutive_noise,
                        repeat_unit, repeat_count, repeat_noise, long_distance_repeats_min_length,
                        long_distance_repeats_noise
                    )
                result_dict[seq_name].extend(
                    [x.add_base(i) for x in _results]
                )
    return result_dict


if __name__ == '__main__':
    # test_seq = 'GTATATGCATGTGATATATATGG'
    # for el in long_distance_segment_repeats_check(test_seq, 5):
    #     print(el.dict())
    # print(test_seq[el.start:el.end])
    test_file = '../tests/data/test1.fasta'
    with open(test_file, 'r') as f:
        for seq_name, results in check_fasta(f.read(), 0.4, 0.6, 15, 4, 0.1, 'AGC', 3, 0.1, 5, 0.1).items():
            df = pd.DataFrame([el.dict() for el in results])
            df.to_excel(f'../tests/data/{seq_name}.xlsx', index=False)
