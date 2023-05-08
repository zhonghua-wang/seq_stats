from collections import Counter
from typing import List, Tuple
import re

from pydantic import BaseModel


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
    seq = seq.upper()
    return (seq.count('G') + seq.count('C')) / len(seq)


# verify the sequence is valid, which contains only 'A', 'T', 'C', 'G'
def is_valid(seq: str) -> bool:
    return re.match(r'^[ATCG]+$', seq, re.IGNORECASE) is not None


# find the regions that GC contents < 40% or > 60% in the specified sequence
def gc_check(seq: str, window_size: int = 15, gc_threshold=(0.4, 0.6)) -> List[CheckResult]:
    seq = seq.upper()
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
        gc_regions.append(CheckResult(start=0, end=end, type='gc', problem=f'GC content < {lower}'))
        gc_regions.extend([x.add_base(end) for x in gc_check(seq[end:], window_size, gc_threshold)])
    elif gc_content_ > upper:
        gc_regions.append(CheckResult(start=0, end=end, type='gc', problem=f'GC content > {upper}'))
        gc_regions.extend([x.add_base(end) for x in gc_check(seq[end:], window_size, gc_threshold)])
    else:
        gc_regions.extend([x.add_base(1) for x in gc_check(seq[1:], window_size, gc_threshold)])
    return gc_regions


# find the regions that contain consecutive nucleotides in the specified sequence
def consecutive_check(seq: str, window_size: int = 4, noise: float = 0.0) -> List[CheckResult]:
    seq = seq.upper()
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
            CheckResult(start=0, end=end, type='consecutive', problem=f'consecutive nucleotides: {mc[0]}')
        )
        consecutive_regions.extend([x.add_base(end) for x in consecutive_check(seq[end:], window_size, noise)])
    else:
        consecutive_regions.extend([x.add_base(1) for x in consecutive_check(seq[1:], window_size, noise)])
    return consecutive_regions


# consecutive repeat unit check
def repeat_unit_check_rex(seq: str, repeat_unit: str, repeat_count: int) -> List[CheckResult]:
    # use regex to find the repeat unit
    pattern = f"(AGC){{{repeat_count},}}"  # 匹配 AGC 重复 3 次以上的模式
    match_list = []
    for match in re.finditer(pattern, seq):
        match_list.append(
            CheckResult(start=match.start(), end=match.end(), type='unit_repeat', problem=f'repeat unit: {repeat_unit}')
        )
    return match_list


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
                type='unit_repeat', problem=f'repeat unit: {repeat_unit}'
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


def long_distance_fixed_segment_repeats_check(dna_sequence, repeat_length):
    """
    在DNA序列中查找多次出现的相同片段序列组以及它们的索引位置。

    参数：
    dna_sequence (str): DNA序列字符串。
    repeat_length (int): 相同片段序列的长度。

    返回：
    repeats (dict): 包含相同片段序列组和索引位置的字典，键为相同片段序列，值为索引位置的列表。
    """
    repeats = {}
    for i in range(len(dna_sequence) - repeat_length + 1):
        repeat = dna_sequence[i:i + repeat_length]
        if repeat in repeats:
            pos_list = repeats[repeat]
            # append to repeats if the repeat is not adjacent to the previous one
            if i - pos_list[-1] - repeat_length > 1:
                repeats[repeat].append(i)
        else:
            repeats[repeat] = [i]

    # 只保留重复出现的片段序列
    result_list = []
    for repeat, indexes in repeats.items():
        if len(indexes) > 1:
            for index in indexes:
                result_list.append(
                    CheckResult(
                        start=index, end=index + repeat_length,
                        type='segment_repeat', problem=f'repeat unit: {repeat}'
                    )
                )
    return result_list


def long_distance_segment_repeats_check(seq: str, min_length: int = 5):
    result_list = []
    if len(seq) < min_length * 2:
        return result_list
    for length in range(int(len(seq) / 2), min_length - 1, -1):
        for result in long_distance_fixed_segment_repeats_check(seq, length):
            if any(
                    [result.start <= x.start <= result.end or result.start <= x.end <= result.end for x in result_list]
            ):
                continue
            result_list.append(result)
    return result_list


if __name__ == '__main__':
    test_seq = 'GTATATGCATGTGATATATATGG'
    for el in long_distance_segment_repeats_check(test_seq, 5):
        print(el)
        # print(test_seq[el.start:el.end])
