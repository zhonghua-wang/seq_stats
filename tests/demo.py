import re

def identify_specific_segments(dna_sequence):
    """
    识别DNA序列中满足条件的特定连续片段

    参数：
    dna_sequence (str): 输入的DNA序列

    返回：
    List[str]: 满足条件的特定连续片段列表
    """
    segments = []
    current_segment = ""
    for base in dna_sequence:
        if base.upper() in ["A", "T", "C", "G"]:
            current_segment += base.upper()
        else:
            if len(current_segment) >= 15:
                gc_content = (current_segment.count("G") + current_segment.count("C")) / len(current_segment)
                if gc_content < 0.4 or gc_content > 0.6:
                    segments.append(current_segment)
            current_segment = ""

    if len(current_segment) >= 15:
        gc_content = (current_segment.count("G") + current_segment.count("C")) / len(current_segment)
        if gc_content < 0.4 or gc_content > 0.6:
            segments.append(current_segment)

    return segments

dna = 'AAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTGGGGGGGGCCCCCCCGGGGGGGGGGGGTC'
gc_regions = identify_specific_segments(dna)
print(gc_regions)