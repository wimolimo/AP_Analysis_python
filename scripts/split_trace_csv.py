"""Convert time,species,... CSV to traces, compositions, and a masslist.

Edit INPUT_CSV below, then press Play in your editor.
Optional command-line use: python scripts/split_trace_csv.py input.csv
Naive input times are interpreted using TIMEZONE below.
Formula masses are calculated from the label as written (no inferred adduct).
"""
import argparse
from collections import Counter
from datetime import date
from pathlib import Path
import re

import pandas as pd


# Settings for running with the editor's Play button.
INPUT_CSV = Path(r"C:\Users\c7441399\Documents\CLOUD18\AP_run\LTOF_UFRA\243K_20251128_20251201_CLOUD18Nonanal_LTOF_UFRA_20260424.csv")
UNITS = "unknown"
TIMEZONE = "UTC"
TRANSMISSION = None  # Unknown; set a number only if known.
OVERWRITE_OUTPUTS = True  # Replace the three generated outputs when rerunning.


MASSES = {"C": 12.0, "C(13)": 13.00335483507, "H": 1.00782503223,
          "N": 14.00307400443, "N(15)": 15.00010889888,
          "O": 15.99491461957, "O(18)": 17.99915961286,
          "S": 31.9720711744, "S(34)": 33.967867004,
          "I": 126.9044719, "Si": 27.97692653465,
          "F": 18.99840316273, "Cr": 51.94050623}
ELEMENTS = ["C", "C(13)", "H", "H+", "N", "O", "O(18)", "S",
            "N(15)", "F", "Cr", "S(34)", "I", "Si"]
MASSLIST_ELEMENTS = ["C", "C(13)", "H", "H+", "N", "O", "O(18)",
                     "S", "S(34)", "I", "Si"]


def composition(label):
    """Parse repeated elements, nested groups, and bracketed isotope labels."""
    text = label.rstrip("+-")
    pos = 0

    def group():
        nonlocal pos
        counts = Counter()
        while pos < len(text) and text[pos] != ")":
            if text[pos] == "(":
                pos += 1
                atoms = group()
                if pos == len(text):
                    raise ValueError(f"Unclosed group: {label}")
                pos += 1
            else:
                match = re.match(r"\[(\d+)([A-Z][a-z]?)\]|([A-Z][a-z]?)", text[pos:])
                if not match:
                    raise ValueError(f"Unsupported formula: {label}")
                isotope, symbol, ordinary = match.groups()
                atoms = {f"{symbol}({isotope})" if isotope else ordinary: 1}
                pos += len(match[0])
            match = re.match(r"\d+", text[pos:])
            multiplier = int(match[0]) if match else 1
            if match:
                pos += len(match[0])
            for atom, count in atoms.items():
                counts[atom] += count * multiplier
        return counts

    try:
        mass = float(label)
        return [float("nan")] * len(ELEMENTS) + [mass, label]
    except ValueError:
        counts = group()
        if pos != len(text) or not counts or set(counts) - MASSES.keys():
            raise ValueError(f"Unsupported formula: {label}")
        mass = sum(MASSES[atom] * count for atom, count in counts.items())
        # Electron mass correction only when a charge is explicit in the label.
        mass += 0.000548579909065 * (label.endswith("-") - label.endswith("+"))
        return [counts.get(atom, 0) for atom in ELEMENTS] + [mass, label]


def write_masslist(path, comp, overwrite=False):
    """Write the example's # Elements / #Masses layout, tab-separated."""
    # Add elements needed by this dataset instead of dropping isotope/F/Cr counts.
    elements = MASSLIST_ELEMENTS + [
        atom for atom in ELEMENTS if atom not in MASSLIST_ELEMENTS
        and comp[atom].fillna(0).ne(0).any()
    ]
    table = comp.reindex(columns=[*elements, "Mass", "SumFormula"], fill_value=0)
    table = table.rename(columns={"SumFormula": "Name"}).copy()
    # Keep unidentified peaks at their supplied mass, matching the example layout.
    unknown = comp[ELEMENTS].isna().all(axis=1)
    table[elements] = table[elements].fillna(0).astype(float)
    table.loc[unknown, "Name"] = ""
    table = table.sort_values("Mass", kind="stable")
    with path.open("w" if overwrite else "x", encoding="utf-8", newline="") as stream:
        stream.write("# Elements: \n" + "\t".join(elements) + "\t\n\n")
        stream.write("\t".join(str(1.007276466621 if e == "H+" else MASSES[e])
                               for e in elements) + "\t\n\n#Masses: \n")
        table.to_csv(stream, sep="\t", index=False)


def write_table(path, table, metadata, overwrite=False):
    lines = [f"{key}:\t{value}" for key, value in metadata.items()]
    with path.open("w" if overwrite else "x", encoding="utf-8", newline="") as stream:
        stream.write(f"number of header rows:\t{len(lines) + 1}\n")
        stream.write("\n".join(lines) + "\n")
        table.to_csv(stream, sep="\t", index=False, na_rep="NaN")


def convert(source, units="unknown", timezone="UTC", transmission=None, overwrite=False):
    source = Path(source)
    # Check the original header before pandas can rename duplicate labels.
    import csv
    with source.open(encoding="utf-8-sig", newline="") as stream:
        labels = next(csv.reader(stream))
    if len(labels) != len(set(labels)):
        raise ValueError("Duplicate column labels: make channel names unique first")
    if labels[0].lower() != "time" or len(labels) < 2:
        raise ValueError("Expected a time column followed by species columns")
    traces = pd.read_csv(source, dtype=str, keep_default_na=False)
    if traces.empty:
        raise ValueError("Input contains no data rows")
    times = pd.to_datetime(traces.iloc[:, 0], errors="raise")
    if times.isna().any():
        raise ValueError("Missing timestamps")
    if times.dt.tz is None:
        times = times.dt.tz_localize(timezone, ambiguous="raise", nonexistent="raise")
    times = times.dt.tz_convert("UTC")
    traces.columns = ["Time", *labels[1:]]
    traces["Time"] = times.dt.strftime("%Y-%m-%dT%H:%M:%S.%f")
    traces.insert(1, "unixTime", [t.timestamp() for t in times])
    comp = pd.DataFrame([composition(s) for s in labels[1:]],
                        columns=[*ELEMENTS, "Mass", "SumFormula"])
    comp["InletTransmission"] = transmission
    trace_path = source.with_name(f"LTOFtraces_{source.stem}.csv")
    comp_path = source.with_name(f"LTOFcompositions_{source.stem}.txt")
    masslist_path = source.with_name(f"masslist_nonanal_LTOF.csv")
    if not overwrite and any(p.exists() for p in (trace_path, comp_path, masslist_path)):
        raise FileExistsError("Output already exists; rename or remove it before rerunning")
    common = {"date": date.today().isoformat(), "level": 1, "version number": "01",
              "column delimiter": "tab"}
    write_table(trace_path, traces, {**common, "title": f"Traces from {source.name}",
                "time format": "Time in UTC ISO 8601; unixTime in seconds",
                "start time": traces.Time[times.idxmin()],
                "end time": traces.Time[times.idxmax()], "column units": units,
                "comment": "Values preserved; no calibration or inlet correction applied."}, overwrite)
    write_table(comp_path, comp, {**common, "title": f"Compositions from {source.name}",
                "column units": "Mass in amu; elemental counts and transmission unitless",
                "comment": "Mass calculated from labels; numeric labels retain their mass with unknown counts. "
                "H+ is zero (no inferred proton adduct). Missing transmission is NaN."}, overwrite)
    write_masslist(masslist_path, comp, overwrite)
    return trace_path, comp_path, masslist_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", nargs="?", type=Path, default=INPUT_CSV)
    parser.add_argument("--units", default=UNITS)
    parser.add_argument("--timezone", default=TIMEZONE)
    parser.add_argument("--transmission", type=float, default=TRANSMISSION,
                        help="unknown transmission; omitted means unknown")
    args = parser.parse_args()
    for output in convert(args.input_csv, args.units, args.timezone, args.transmission,
                          overwrite=OVERWRITE_OUTPUTS):
        print(output)
