from seq_stats.check import *


def test_gc_check():
    seq1 = 'TGCTACATTGTGCTGCGTTT'
    seq2 = 'TAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTGGGGGGGGCCCCCCCGGGGGGGGGGGGTCGGGGGGGCCCCCCTGCT'
    assert gc_check(seq1) == []
    result2 = gc_check(seq2)
    assert len(result2) == 2
    assert result2[0].start == 0 and result2[0].end == 63
    assert result2[0].type == 'gc' and result2[0].problem == 'GC content < 0.4'


def test_consecutive_check():
    seq1 = 'TGCTACATTGTGCTGCGTTT'
    seq2 = 'TAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTGGGGGGGGCCCCCCCGGGGGGGGGGGGTCGGGGGGGCCCCCCTGCT'
    seq3 = 'ATTTGTTTTCCCTTGCG'
    assert consecutive_check(seq1) == []
    result2 = consecutive_check(seq2)
    assert len(result2) == 7
    result3 = consecutive_check(seq3)
    result4 = consecutive_check(seq3, 4, 1 / 7)
    assert len(result3) == 1
    assert len(result4) == 1
    assert result3[0].start == 5 and result3[0].end == 9
    assert result4[0].start == 1 and result4[0].end == 9
    assert seq3[result4[0].start:result4[0].end] == 'TTTGTTTT'


def test_repeat_unit_check_rex():
    seq1 = 'AGAAGCAGCTTA'
    seq2 = 'AGAGCAGCAGCAGCTAGCAGCAGCTTTAAA'
    result1 = repeat_unit_check_rex(seq1, 'AGC', 3)
    result2 = repeat_unit_check_rex(seq2, 'AGC', 3)
    assert len(result1) == 0
    assert len(result2) == 2
    assert seq2[result2[0].start: result2[0].end] == 'AGC' * 4


def test_repeat_unit_check():
    seq1 = 'AGAAGCAGCTTA'
    seq2 = 'AGAGCAGCAGCAGCTAGCAGCAGCTTTAAA'
    result1 = repeat_unit_check(seq1, 'AGC', 3)
    result2 = repeat_unit_check(seq2, 'AGC', 3)
    result3 = repeat_unit_check(seq2, 'AGC', 3, 1 / 7)
    result4 = repeat_unit_check(seq2, 'AGC', 3, 1 / 8)
    assert len(result1) == 0
    assert len(result2) == 2
    assert seq2[result2[0].start: result2[0].end] == 'AGC' * 4
    assert len(result3) == 1
    assert result3[0].start == 2 and result3[0].end == 24
    assert len(result4) == 2


def test_long_distance_segment_repeats_check():
    test_seq = 'GTATATGCATGTGATATATATGG'
    result = long_distance_segment_repeats_check(test_seq, 5)
    assert len(result) == 2
