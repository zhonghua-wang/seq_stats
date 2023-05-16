import os
import sys
from time import sleep

import pandas as pd
from gooey import Gooey, GooeyParser

from seq_stats.check import check_fasta


@Gooey(
    language='chinese', auto_start=True, program_name='序列检查',
    advanced=True,
    progress_regex=r"^进度: (?P<current>\d+)/(?P<total>\d+)$",
    progress_expr="current / total * 100",
    timing_options={
        'show_time_remaining': True,
        'hide_time_remaining_on_complete': True
    }
)
def main():
    parser = GooeyParser()
    parser.add_argument('input', help='输入待分析Fasta文件', widget='FileChooser')
    parser.add_argument('output', help='结果目录', widget='DirChooser')

    gc = parser.add_argument_group('GC含量检查', gooey_options={'show_border': True})
    gc.add_argument('--gc_min', help='最小GC含量', default=0.4, type=float)
    gc.add_argument('--gc_max', help='最大GC含量', default=0.6, type=float)

    consecutive = parser.add_argument_group('单一重复序列检查', gooey_options={'show_border': True})
    consecutive.add_argument('--consecutive_window_size', help='重复序列窗口大小', default=4, type=int)
    consecutive.add_argument('--consecutive_noise', help='重复序列允许噪声比例', default=0, type=float)

    repeat_unit = parser.add_argument_group('连续重复单元检查', gooey_options={'show_border': True})
    repeat_unit.add_argument('--repeat_unit', help='重复单元', default='AGC')
    repeat_unit.add_argument('--repeat_count', help='重复次数', default=3, type=int)
    repeat_unit.add_argument('--repeat_noise', help='重复单元允许噪声比例', default=0, type=float)

    long_distance_repeats = parser.add_argument_group('长距离重复序列检查', gooey_options={'show_border': True})
    long_distance_repeats.add_argument('--long_distance_repeats_min_length', help='最小重复序列长度', default=5,
                                       type=int)
    long_distance_repeats.add_argument('--long_distance_repeats_noise', help='重复序列允许噪声比例', default=0,
                                       type=float)

    args = parser.parse_args()
    file_name = args.input.split('/')[-1].split('.')[0]
    result_file = os.path.join(args.output, f'{file_name}-result.xlsx')

    excel_writer = pd.ExcelWriter(result_file)

    with open(args.input, 'r') as fasta_file:
        for seq_name, results in check_fasta(
            fasta_file.read(),
            args.gc_min, args.gc_max, 15,
            args.consecutive_window_size, args.consecutive_noise,
            args.repeat_unit, args.repeat_count, args.repeat_noise,
            args.long_distance_repeats_min_length, args.long_distance_repeats_noise
        ).items():
            df = pd.DataFrame([el.dict() for el in results])
            df.to_excel(excel_writer, sheet_name=seq_name, index=False)

    excel_writer.save()
    print('分析完成, 结果文件已保存至: ', result_file)


if __name__ == '__main__':
    main()
