"""Remove duplicate/near-duplicate Mass rows, preserving the original text.

Edit INPUT_CSV and MASS_TOLERANCE below, then run this script.
Keeps the first row in file order; compares each later mass with retained rows.
Exactly equal masses are always duplicates. Other differences must be strictly
less than the tolerance (0 means exact duplicates only).
"""
from bisect import bisect_left, insort
import csv
from datetime import datetime
from decimal import Decimal, InvalidOperation
from pathlib import Path
import os
import shutil
import tempfile

# INPUT_CSV = Path(r"C:\Users\c7441399\Documents\CLOUD18\AP_run\masslist_v4.csv")
INPUT_CSV = Path(r"C:\Users\c7441399\Documents\CLOUD18\AP_run\LTOF_UFRA\masslist_nonanal_LTOF.csv")
MASS_TOLERANCE = 0.0001


def remove_duplicate_masses(input_csv, tolerance=0.01):
    """Overwrite input after validation and backup; return removed-row count.

    Supports comma-, tab-, or semicolon-separated tables with optional preamble.
    Invalid/missing masses abort without modifying the input.
    """
    path = Path(input_csv)
    threshold = Decimal(str(tolerance))
    if not threshold.is_finite() or threshold < 0:
        raise ValueError("Tolerance must be finite and nonnegative")
    with path.open(encoding="utf-8-sig", newline="") as stream:
        text = stream.read()
    lines = text.splitlines(keepends=True)
    for header_index, line in enumerate(lines):
        for separator in ("\t", ",", ";"):
            columns = next(csv.reader([line], delimiter=separator))
            if "Mass" in [c.strip() for c in columns]:
                mass_index = [c.strip() for c in columns].index("Mass")
                break
        else:
            continue
        break
    else:
        raise ValueError("No table header containing a Mass column found")

    kept = lines[:header_index + 1]
    masses = []
    retained_rows = {}
    deletion_details = []
    removed = 0
    # reader.line_num also preserves quoted records spanning several lines.
    reader = csv.reader(lines[header_index + 1:], delimiter=separator)
    previous_line = 0
    for row in reader:
        line_number = header_index + 2 + previous_line
        raw = lines[header_index + 1 + previous_line:header_index + 1 + reader.line_num]
        previous_line = reader.line_num
        if not row or not any(cell.strip() for cell in row):
            kept.extend(raw)
            continue
        try:
            mass = Decimal(row[mass_index].strip())
            if not mass.is_finite():
                raise InvalidOperation
        except (InvalidOperation, IndexError):
            raise ValueError(f"Invalid Mass near line {header_index + 1 + reader.line_num}") from None
        index = bisect_left(masses, mass)
        neighbors = masses[max(0, index-1):index+1]
        matches = [other for other in neighbors
                   if mass == other or abs(mass-other) < threshold]
        if matches:
            matched_mass = min(matches, key=lambda other: abs(mass-other))
            deletion_details.append((line_number, "".join(raw).rstrip("\r\n"),
                                     mass, matched_mass, retained_rows[matched_mass]))
            removed += 1
        else:
            insort(masses, mass)
            retained_rows[mass] = (line_number, "".join(raw).rstrip("\r\n"))
            kept.extend(raw)

    if not removed:
        print(f"No duplicate masses found: {path}")
        return 0
    backup = path.with_name(path.name + datetime.now().strftime(".%Y%m%d_%H%M%S_%f.bak"))
    shutil.copy2(path, backup)
    encoding = "utf-8-sig" if path.read_bytes().startswith(b"\xef\xbb\xbf") else "utf-8"
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding=encoding, newline="",
                                         dir=path.parent, delete=False) as stream:
            temporary = stream.name
            stream.write("".join(kept))
        os.replace(temporary, path)
    finally:
        if temporary and os.path.exists(temporary):
            os.unlink(temporary)
    print("\nColumns: " + lines[header_index].rstrip("\r\n"))
    for line_number, deleted_row, mass, matched_mass, retained in deletion_details:
        kept_line, kept_row = retained
        print(f"\nDeleted original line {line_number}: Mass={mass}")
        print(deleted_row)
        print(f"Matched retained original line {kept_line}: Mass={matched_mass}")
        print(kept_row)
        print(f"Absolute mass difference: {abs(mass - matched_mass)}")
    print(f"\nRemoved {removed} duplicate rows; kept {len(masses)} mass rows. Tolerance: {threshold}")
    print(f"Saved: {path}\nBackup: {backup}")
    return removed


if __name__ == "__main__":
    remove_duplicate_masses(INPUT_CSV, MASS_TOLERANCE)
