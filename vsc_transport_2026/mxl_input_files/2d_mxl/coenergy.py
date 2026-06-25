import argparse
import re
from pathlib import Path

import MDAnalysis as mda
import numpy as np


GRID_SHAPE = (144, 144)
N_BOXES = int(np.prod(GRID_SHAPE))
N_STEPS = 401
OUTPUT_NAME = "coenergy.npy"
MOLECULAR_ID_RE = re.compile(r"\[MaxwellLink\] Assigned a molecular ID:\s*(\d+)")


def read_coenergy_file(path):
    """Read the second numeric column from one coenergy_box_*.xyz file."""
    values = []

    with mda.lib.util.openany(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"{path}:{line_number} does not have two columns")

            try:
                values.append(float(parts[1]))
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number} second column is not a number: {parts[1]!r}"
                ) from exc

    values = np.asarray(values, dtype=np.float64)
    if values.shape != (N_STEPS,):
        raise ValueError(
            f"{path} has {values.size} data rows; expected {N_STEPS}"
        )

    return values


def read_molecular_id(path):
    """Read the MaxwellLink molecular/grid ID from one out_*.lmp file."""
    with mda.lib.util.openany(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            match = MOLECULAR_ID_RE.search(line)
            if match:
                return int(match.group(1))

    raise ValueError(
        f"{path} does not contain a '[MaxwellLink] Assigned a molecular ID' line"
    )


def collect_molecular_ids(folder):
    """Map each coenergy_box_N.xyz file index to its assigned flat grid ID."""
    box_to_grid = {}
    grid_to_box = {}
    unresolved_boxes = []

    for box_index in range(1, N_BOXES + 1):
        path = folder / f"out_{box_index}.lmp"
        if not path.is_file():
            unresolved_boxes.append((box_index, f"missing required file: {path}"))
            continue

        try:
            grid_id = read_molecular_id(path)
        except ValueError as exc:
            unresolved_boxes.append((box_index, str(exc)))
            continue

        if not 0 <= grid_id < N_BOXES:
            raise ValueError(
                f"{path} assigns molecular ID {grid_id}; expected 0..{N_BOXES - 1}"
            )
        if grid_id in grid_to_box:
            raise ValueError(
                f"Duplicate molecular ID {grid_id}: "
                f"out_{grid_to_box[grid_id]}.lmp and out_{box_index}.lmp"
            )

        box_to_grid[box_index] = grid_id
        grid_to_box[grid_id] = box_index

    missing_grid_ids = sorted(set(range(N_BOXES)) - set(grid_to_box))

    if unresolved_boxes:
        if len(unresolved_boxes) == 1 and len(missing_grid_ids) == 1:
            box_index, reason = unresolved_boxes[0]
            grid_id = missing_grid_ids[0]
            box_to_grid[box_index] = grid_id
            grid_to_box[grid_id] = box_index
            print(
                "Warning: inferred molecular ID "
                f"{grid_id} for out_{box_index}.lmp because it was the only "
                f"unresolved output file ({reason})."
            )
            missing_grid_ids = []
        else:
            details = "; ".join(
                f"out_{box_index}.lmp: {reason}"
                for box_index, reason in unresolved_boxes
            )
            raise ValueError(
                "Could not infer molecular IDs because "
                f"{len(unresolved_boxes)} output files were unresolved: {details}"
            )

    if missing_grid_ids:
        raise ValueError(f"Missing molecular IDs: {missing_grid_ids}")

    return box_to_grid


def collect_coenergy(folder):
    folder = Path(folder).expanduser().resolve()
    if not folder.is_dir():
        raise NotADirectoryError(f"{folder} is not a directory")

    box_to_grid = collect_molecular_ids(folder)
    data = np.empty((N_BOXES, N_STEPS), dtype=np.float64)
    for box_index in range(1, N_BOXES + 1):
        path = folder / f"coenergy_box_{box_index}.xyz"
        if not path.is_file():
            raise FileNotFoundError(f"Missing required file: {path}")

        data[box_to_grid[box_index]] = read_coenergy_file(path)

    # Molecular IDs are flattened 0-based grid IDs.  Reshaping in C order
    # places ID ``i * GRID_SHAPE[1] + j`` at grid position ``(i, j)``.
    return data.T.reshape((N_STEPS, *GRID_SHAPE))


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Collect coenergy_box_1.xyz ... coenergy_box_20736.xyz into a "
            "1001 x 144 x 144 NumPy array, using out_*.lmp MaxwellLink "
            "molecular IDs to place each file on the correct grid point."
        )
    )
    parser.add_argument(
        "folder",
        help="Folder containing coenergy_box_*.xyz and out_*.lmp files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=OUTPUT_NAME,
        help=f"Output .npy filename, saved inside folder. Default: {OUTPUT_NAME}",
    )
    args = parser.parse_args()

    folder = Path(args.folder).expanduser().resolve()
    output_path = folder / args.output
    data = collect_coenergy(folder)
    np.save(output_path, data)
    print(f"Saved {data.shape} array to {output_path}")


if __name__ == "__main__":
    main()
