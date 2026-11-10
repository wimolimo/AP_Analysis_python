from __future__ import annotations
from pathlib import Path
from datetime import datetime
import re

_STAGE_DATA = []


def parse_cloud_log(filepath):
    """Parse stages, keeping event text and untimed notes in source order.

    Only event timestamps position markers. Untimed text is explicitly labeled
    and never inherits a preceding or following event's time. Metadata fields
    are read on their own lines; decorative headers and staff lists are skipped.
    """
    content = Path(filepath).read_text(encoding="utf-8", errors="replace")
    stage_numbers = re.findall(r"\*{12,}(\d+\.\d+)\*{12,}", content)
    blocks = re.split(r"\*{12,}\d+\.\d+\*{12,}", content)

    results = []
    for idx, number in enumerate(stage_numbers):
        block = blocks[idx + 1] if idx + 1 < len(blocks) else ""
        fields = {"Type": "", "Description": "", "Comments": ""}
        details = []
        stage_time = None
        in_events = False
        for raw_line in block.splitlines():
            line = raw_line.strip()
            if not line:
                continue
            field = re.fullmatch(r"(Type|Description|Comments):[ \t]*(.*)", line)
            if field and not in_events:
                fields[field[1]] = field[2]
                continue
            # Decorative headings and staff information are not event notes.
            if "CLOUD Run.Stage=" in line:
                continue
            if set(line) == {"*"}:
                in_events = True
                continue
            if line.startswith(("Shift team leader:", "Team member:")):
                continue
            timestamp = re.match(r"(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})(.*)", line)
            if timestamp:
                event_time = datetime.strptime(timestamp[1], "%Y-%m-%d %H:%M:%S")
                if stage_time is None:
                    stage_time = event_time
                in_events = True
                text = timestamp[2].strip()
                details.append(timestamp[1] + (" " + text if text else ""))
            elif in_events:
                # Preserve source order without borrowing a neighboring time.
                # A standalone timestamp followed by text remains two entries.
                details.append("[untimed note] " + line)
        results.append({
            "stage": number,
            "time": stage_time,
            "type": fields["Type"],
            "description": fields["Description"],
            "comments": fields["Comments"],
            "details": details,
        })
    return results


def load_stages(filepath):
    """Parse a CLOUD log and store its stages for subsequent plotting."""
    global _STAGE_DATA
    from .prepared import PreparedSource, read_prepared_stages
    _STAGE_DATA = (read_prepared_stages(filepath) if isinstance(filepath, PreparedSource)
                   else parse_cloud_log(filepath))
    print(f"Loaded {len(_STAGE_DATA)} stages")
    return _STAGE_DATA


def get_stages():
    """Return the stages most recently loaded with :func:`load_stages`."""
    return _STAGE_DATA
