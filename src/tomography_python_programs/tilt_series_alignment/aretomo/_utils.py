from typing import Tuple

def coerce_gpu_ids(
        gpu_ids: str
    ) -> Tuple:
    """Convert GPU  ID input from colon spaced (1:2:3) to no space inbetween GPU IDs (123)"""
    return tuple(map(int, gpu_ids.strip().replace(':','')))
