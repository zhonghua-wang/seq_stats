import tempfile
from typing import List, Annotated

from fastapi import FastAPI, UploadFile
from fast_ss import settings
from seq_stats.check import *

app = FastAPI()


@app.post("/check_seq")
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


@app.post("/check_file/")
async def check_file(
        file: UploadFile,
        gc_min: float, gc_max: float, gc_window_size: int,
        consecutive_window_size: int, consecutive_noise: float,
        repeat_unit: str, repeat_count: int, repeat_noise: float,
        long_distance_repeats_min_length: int, long_distance_repeats_noise: float
):
    contents = await file.read()
    fa = contents.decode('utf-8').split('\n')
    seq = ''.join(fa[1:])
    return check_seq(
        seq,
        gc_min, gc_max, gc_window_size,
        consecutive_window_size, consecutive_noise,
        repeat_unit, repeat_count, repeat_noise,
        long_distance_repeats_min_length, long_distance_repeats_noise
    )


@app.post("/check_files/")
async def check_files(
        files: list[UploadFile],
        # gc_min: float, gc_max: float, gc_window_size: int,
        # consecutive_window_size: int, consecutive_noise: float,
        # repeat_unit: str, repeat_count: int, repeat_noise: float,
        # long_distance_repeats_min_length: int, long_distance_repeats_noise: float
):
    print(settings.Settings().media_path)
    print(files)
    return
        # create temp dir and package files to a zip file

