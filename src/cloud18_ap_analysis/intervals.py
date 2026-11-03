"""Shared UTC averaging intervals. Edit here to update all three run scripts.

List order defines t_i; keep excluded intervals commented out.
"""
from datetime import datetime

important_time_ranges_m50C = [
    (datetime(2025, 9, 25, 10, 0, 0), datetime(2025, 9, 25, 20, 0, 0)),
    (datetime(2025, 9, 26, 1, 30, 0), datetime(2025, 9, 26, 6, 0, 0)),
    (datetime(2025, 9, 26, 12, 0, 0), datetime(2025, 9, 26, 16, 0, 0)),
    (datetime(2025, 9, 26, 19, 0, 0), datetime(2025, 9, 27, 0, 0, 0)),
    (datetime(2025, 9, 27, 4, 0, 0), datetime(2025, 9, 27, 4, 30, 0)),
]

important_time_ranges_m25C = [
    (datetime(2025, 9, 30, 2, 40, 0), datetime(2025, 9, 30, 3, 40, 0)),
    (datetime(2025, 9, 30, 8, 10, 0), datetime(2025, 9, 30, 9, 10, 0)),
    # (datetime(2025, 9, 30, 11, 40, 0), datetime(2025, 9, 30, 12, 30, 0)),     # fan speed at 100%
    (datetime(2025, 9, 30, 14, 30, 0), datetime(2025, 9, 30, 15, 20, 0)),
    (datetime(2025, 9, 30, 16, 10, 0), datetime(2025, 9, 30, 17, 0, 0)),
    (datetime(2025, 9, 30, 18, 30, 0), datetime(2025, 9, 30, 21, 0, 0)),
    (datetime(2025, 9, 30, 22, 0, 0), datetime(2025, 9, 30, 23, 30, 0)),
    (datetime(2025, 10, 1, 2, 10, 0), datetime(2025, 10, 1, 3, 15, 0)),
    (datetime(2025, 10, 1, 6, 0, 0), datetime(2025, 10, 1, 7, 30, 0)),
    (datetime(2025, 10, 1, 15, 20, 0), datetime(2025, 10, 1, 16, 20, 0)),
    (datetime(2025, 10, 1, 21, 30, 0), datetime(2025, 10, 1, 22, 15, 0)),
    (datetime(2025, 10, 2, 3, 0, 0), datetime(2025, 10, 2, 4, 15, 0)),
    (datetime(2025, 10, 2, 5, 10, 0), datetime(2025, 10, 2, 5, 30, 0)),
    (datetime(2025, 10, 2, 7, 45, 0), datetime(2025, 10, 2, 9, 0, 0)),
    (datetime(2025, 10, 2, 11, 7, 0), datetime(2025, 10, 2, 11, 35, 0)),
]

important_time_ranges_m10C = [
    # (datetime(2025, 10, 8, 18, 0, 0), datetime(2025, 10, 8, 19, 0, 0)),      
    # (datetime(2025, 10, 8, 20, 50, 0), datetime(2025, 10, 8, 21, 50, 0)),    
    (datetime(2025, 10, 9, 2, 30, 0), datetime(2025, 10, 9, 3, 0, 0)),
    (datetime(2025, 10, 9, 3, 30, 0), datetime(2025, 10, 9, 4, 15, 0)),
    (datetime(2025, 10, 9, 5, 0, 0), datetime(2025, 10, 9, 5, 30, 0)),          # split this stage to due O3 variations
    (datetime(2025, 10, 9, 6, 0, 0), datetime(2025, 10, 9, 6, 30, 0)),
    (datetime(2025, 10, 9, 6, 50, 0), datetime(2025, 10, 9, 7, 15, 0)),
    (datetime(2025, 10, 9, 7, 30, 0), datetime(2025, 10, 9, 8, 0, 0)),
    (datetime(2025, 10, 9, 15, 30, 0), datetime(2025, 10, 9, 16, 5, 0)),
    (datetime(2025, 10, 9, 23, 15, 0), datetime(2025, 10, 10, 0, 15, 0)),
    (datetime(2025, 10, 10, 1, 15, 0), datetime(2025, 10, 10, 2, 15, 0)),
    (datetime(2025, 10, 10, 5, 30, 0), datetime(2025, 10, 10, 6, 10, 0)),
    (datetime(2025, 10, 10, 6, 35, 0), datetime(2025, 10, 10, 7, 30, 0)),
    (datetime(2025, 10, 10, 9, 30, 0), datetime(2025, 10, 10, 10, 35, 0)),
    (datetime(2025, 10, 10, 10, 50, 0), datetime(2025, 10, 10, 11, 40, 0)),
    (datetime(2025, 10, 10, 12, 50, 0), datetime(2025, 10, 10, 13, 20, 0)),
    (datetime(2025, 10, 10, 15, 0, 0), datetime(2025, 10, 10, 16, 0, 0)),
    (datetime(2025, 10, 10, 16, 50, 0), datetime(2025, 10, 10, 17, 40, 0)),     # due to NH3 changes
    (datetime(2025, 10, 10, 19, 0, 0), datetime(2025, 10, 10, 19, 30, 0)),
    (datetime(2025, 10, 10, 19, 50, 0), datetime(2025, 10, 10, 20, 10, 0)),
]

important_time_ranges_m50C_2 = [
    (datetime(2025, 10, 12, 9, 0, 0), datetime(2025, 10, 12, 10, 5, 0)),
    (datetime(2025, 10, 12, 11, 15, 0), datetime(2025, 10, 12, 12, 10, 0)),
    (datetime(2025, 10, 12, 13, 50, 0), datetime(2025, 10, 12, 14, 30, 0)),
    (datetime(2025, 10, 12, 14, 50, 0), datetime(2025, 10, 12, 15, 8, 0)),
    (datetime(2025, 10, 13, 13, 30, 0), datetime(2025, 10, 13, 14, 50, 0)),
    (datetime(2025, 10, 13, 16, 0, 0), datetime(2025, 10, 13, 16, 30, 0)),
]

IMPORTANT_TIME_RANGES = {
    "m50_1": important_time_ranges_m50C,
    "m25": important_time_ranges_m25C,
    "m10": important_time_ranges_m10C,
    "m50_2": important_time_ranges_m50C_2,
}
