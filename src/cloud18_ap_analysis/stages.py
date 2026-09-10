from __future__ import annotations
from pathlib import Path
from datetime import datetime
import re

_STAGE_DATA = []


def parse_cloud_log(filepath):
    content = Path(filepath).read_text(encoding="utf-8", errors="replace")
    stage_numbers = re.findall(r"\*{12,}(\d+\.\d+)\*{12,}", content)
    blocks = re.split(r"\*{12,}\d+\.\d+\*{12,}", content)

    results = []
    for idx, number in enumerate(stage_numbers):
        block = blocks[idx + 1] if idx + 1 < len(blocks) else ""
        type_m = re.search(r"Type:\s*(.+)", block)
        desc_m = re.search(r"Description:\s*(.*)", block)
        comm_m = re.search(r"Comments:\s*(.*)", block)
        details = [m.group(0).strip() for m in re.finditer(
            r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} .+", block
        )]
        time_m = re.search(r"(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})", block)
        stage_time = (
            datetime.strptime(time_m.group(1), "%Y-%m-%d %H:%M:%S")
            if time_m else None
        )
        results.append({
            "stage": number,
            "time": stage_time,
            "type": type_m.group(1).strip() if type_m else "",
            "description": desc_m.group(1).strip() if desc_m else "",
            "comments": comm_m.group(1).strip() if comm_m else "",
            "details": details,
        })
    return results


def load_stages(filepath):
    global _STAGE_DATA
    _STAGE_DATA = parse_cloud_log(filepath)
    print(f"Loaded {len(_STAGE_DATA)} stages")
    return _STAGE_DATA


def get_stages():
    return _STAGE_DATA
