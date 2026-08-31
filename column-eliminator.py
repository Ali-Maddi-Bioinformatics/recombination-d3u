#!/usr/bin/env python3

import argparse
import csv
import os
import sys
from pathlib import Path
from typing import List, Tuple

from tqdm import tqdm


def parse_column_numbers(values: List[str]) -> List[int]:
    """
    Parse column numbers.

    Supports:
        2 5 7
        2,5,7
        2,5 7

    Column numbers are 1-based.
    """
    columns = []

    for value in values:
        parts = value.split(",")

        for part in parts:
            part = part.strip()

            if not part:
                continue

            try:
                column_number = int(part)
            except ValueError:
                raise ValueError(
                    f"Invalid column number: {part!r}. "
                    "Column numbers must be integers."
                )

            if column_number < 1:
                raise ValueError(
                    f"Invalid column number: {column_number}. "
                    "Column numbering starts from 1."
                )

            columns.append(column_number)

    # Remove duplicates while preserving numerical order.
    return sorted(set(columns))


def build_output_path(input_path: Path) -> Path:
    """
    Example:
        data.csv -> data_final.csv
    """
    return input_path.with_name(
        f"{input_path.stem}_final{input_path.suffix}"
    )


def build_temp_path(output_path: Path) -> Path:
    """
    Temporary file used until verification succeeds.
    """
    return output_path.with_name(
        f".{output_path.stem}.tmp{output_path.suffix}"
    )


def get_csv_info(
    file_path: Path,
    encoding: str
) -> Tuple[List[str], int]:
    """
    Read CSV header and count data rows.

    Returns:
        header
        data_row_count
    """
    with file_path.open(
        "r",
        encoding=encoding,
        newline=""
    ) as f:

        reader = csv.reader(f)

        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("CSV file is empty.")

        row_count = 0

        print("\n[1/3] Scanning input CSV...")

        for _ in tqdm(
            reader,
            desc="Counting rows",
            unit="rows",
            dynamic_ncols=True
        ):
            row_count += 1

    return header, row_count


def validate_columns(
    header: List[str],
    columns_to_remove: List[int]
) -> List[int]:
    """
    Validate requested 1-based column numbers
    and convert them to zero-based indexes.
    """
    column_count = len(header)

    invalid_columns = [
        col
        for col in columns_to_remove
        if col > column_count
    ]

    if invalid_columns:
        raise ValueError(
            f"Invalid column number(s): {invalid_columns}\n"
            f"The CSV contains only {column_count} columns."
        )

    if len(columns_to_remove) == column_count:
        raise ValueError(
            "You requested removal of all columns. "
            "At least one column must remain."
        )

    return [col - 1 for col in columns_to_remove]


def create_filtered_csv(
    input_path: Path,
    temp_output_path: Path,
    indexes_to_remove: List[int],
    total_rows: int,
    encoding: str
):
    """
    Create a new CSV without selected columns.
    """
    remove_set = set(indexes_to_remove)

    print("\n[2/3] Creating filtered CSV...")

    with input_path.open(
        "r",
        encoding=encoding,
        newline=""
    ) as input_file, temp_output_path.open(
        "w",
        encoding=encoding,
        newline=""
    ) as output_file:

        reader = csv.reader(input_file)
        writer = csv.writer(output_file)

        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("Input CSV is empty.")

        filtered_header = [
            value
            for index, value in enumerate(header)
            if index not in remove_set
        ]

        writer.writerow(filtered_header)

        for row_number, row in enumerate(
            tqdm(
                reader,
                total=total_rows,
                desc="Removing columns",
                unit="rows",
                dynamic_ncols=True
            ),
            start=2
        ):
            if len(row) != len(header):
                raise ValueError(
                    f"Malformed CSV row at line/record {row_number}: "
                    f"expected {len(header)} columns, "
                    f"found {len(row)}."
                )

            filtered_row = [
                value
                for index, value in enumerate(row)
                if index not in remove_set
            ]

            writer.writerow(filtered_row)


def verify_files(
    input_path: Path,
    output_path: Path,
    indexes_to_remove: List[int],
    total_rows: int,
    encoding: str
):
    """
    Verify that:

    1. Only requested columns were removed.
    2. Header is correct.
    3. No data rows were removed.
    4. No data rows were added.
    5. Row order is unchanged.
    6. Remaining field values are exactly identical.
    7. Remaining column order is unchanged.
    """
    remove_set = set(indexes_to_remove)

    print("\n[3/3] Verifying output CSV...")

    with input_path.open(
        "r",
        encoding=encoding,
        newline=""
    ) as original_file, output_path.open(
        "r",
        encoding=encoding,
        newline=""
    ) as generated_file:

        original_reader = csv.reader(original_file)
        generated_reader = csv.reader(generated_file)

        try:
            original_header = next(original_reader)
        except StopIteration:
            raise RuntimeError("Original CSV unexpectedly became empty.")

        try:
            generated_header = next(generated_reader)
        except StopIteration:
            raise RuntimeError("Generated CSV is empty.")

        expected_header = [
            value
            for index, value in enumerate(original_header)
            if index not in remove_set
        ]

        if generated_header != expected_header:
            raise RuntimeError(
                "Verification failed: output header does not match "
                "the expected header."
            )

        original_row_count = 0
        generated_row_count = 0

        with tqdm(
            total=total_rows,
            desc="Verifying rows",
            unit="rows",
            dynamic_ncols=True
        ) as progress:

            while True:
                try:
                    original_row = next(original_reader)
                    original_exists = True
                except StopIteration:
                    original_row = None
                    original_exists = False

                try:
                    generated_row = next(generated_reader)
                    generated_exists = True
                except StopIteration:
                    generated_row = None
                    generated_exists = False

                # Both files ended at the same point.
                if not original_exists and not generated_exists:
                    break

                if original_exists:
                    original_row_count += 1

                if generated_exists:
                    generated_row_count += 1

                # One file ended before the other.
                if original_exists != generated_exists:
                    raise RuntimeError(
                        "Verification failed: row count is different "
                        "between original and generated CSV."
                    )

                if len(original_row) != len(original_header):
                    raise RuntimeError(
                        f"Verification failed: original row "
                        f"{original_row_count + 1} contains "
                        f"{len(original_row)} columns instead of "
                        f"{len(original_header)}."
                    )

                expected_row = [
                    value
                    for index, value in enumerate(original_row)
                    if index not in remove_set
                ]

                if generated_row != expected_row:
                    raise RuntimeError(
                        "\nVerification failed!\n"
                        f"Mismatch found at data row "
                        f"{original_row_count}.\n\n"
                        f"Expected:\n{expected_row}\n\n"
                        f"Generated:\n{generated_row}"
                    )

                progress.update(1)

        if original_row_count != total_rows:
            raise RuntimeError(
                "Verification failed: original row count changed "
                "during processing."
            )

        if generated_row_count != total_rows:
            raise RuntimeError(
                f"Verification failed: expected {total_rows} rows "
                f"but generated file contains "
                f"{generated_row_count} rows."
            )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Remove selected columns from a CSV file and verify "
            "the generated file against the original."
        )
    )

    parser.add_argument(
        "csv_file",
        help="Path to the input CSV file."
    )

    parser.add_argument(
        "columns",
        nargs="+",
        help=(
            "1-based column numbers to remove. "
            "Examples: 2 5 7   or   2,5,7"
        )
    )

    parser.add_argument(
        "--encoding",
        default="utf-8-sig",
        help=(
            "CSV encoding. Default: utf-8-sig "
            "(compatible with UTF-8 and Excel UTF-8 CSV files)."
        )
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite an existing *_final.csv file."
    )

    args = parser.parse_args()

    input_path = Path(args.csv_file).expanduser().resolve()

    if not input_path.exists():
        print(
            f"\nERROR: Input file does not exist:\n{input_path}",
            file=sys.stderr
        )
        sys.exit(1)

    if not input_path.is_file():
        print(
            f"\nERROR: The provided path is not a file:\n{input_path}",
            file=sys.stderr
        )
        sys.exit(1)

    if input_path.suffix.lower() != ".csv":
        print(
            f"\nERROR: Input file must have a .csv extension:\n"
            f"{input_path}",
            file=sys.stderr
        )
        sys.exit(1)

    try:
        columns_to_remove = parse_column_numbers(args.columns)
    except ValueError as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)

    if not columns_to_remove:
        print(
            "\nERROR: No columns were specified.",
            file=sys.stderr
        )
        sys.exit(1)

    output_path = build_output_path(input_path)
    temp_output_path = build_temp_path(output_path)

    if output_path == input_path:
        print(
            "\nERROR: Output path cannot be the same as input path.",
            file=sys.stderr
        )
        sys.exit(1)

    if output_path.exists() and not args.overwrite:
        print(
            f"\nERROR: Output file already exists:\n"
            f"{output_path}\n\n"
            f"Use --overwrite if you want to replace it.",
            file=sys.stderr
        )
        sys.exit(1)

    # Remove stale temporary file if present.
    if temp_output_path.exists():
        temp_output_path.unlink()

    try:
        header, total_rows = get_csv_info(
            input_path,
            args.encoding
        )

        indexes_to_remove = validate_columns(
            header,
            columns_to_remove
        )

        print("\nCSV information")
        print("-" * 60)
        print(f"Input file       : {input_path}")
        print(f"Total columns    : {len(header):,}")
        print(f"Total data rows  : {total_rows:,}")
        print(
            f"Columns to remove: "
            f"{', '.join(map(str, columns_to_remove))}"
        )

        print("\nSelected columns:")
        for column_number in columns_to_remove:
            print(
                f"  {column_number:>5} -> "
                f"{header[column_number - 1]}"
            )

        print(
            f"\nRemaining columns: "
            f"{len(header) - len(columns_to_remove):,}"
        )

        create_filtered_csv(
            input_path=input_path,
            temp_output_path=temp_output_path,
            indexes_to_remove=indexes_to_remove,
            total_rows=total_rows,
            encoding=args.encoding
        )

        verify_files(
            input_path=input_path,
            output_path=temp_output_path,
            indexes_to_remove=indexes_to_remove,
            total_rows=total_rows,
            encoding=args.encoding
        )

        # Verification succeeded.
        # Only now replace/create the final output file.
        os.replace(temp_output_path, output_path)

        print("\n" + "=" * 60)
        print("SUCCESS")
        print("=" * 60)
        print("CSV processing and verification completed successfully.")
        print()
        print(f"Original file : {input_path}")
        print(f"Output file   : {output_path}")
        print(f"Rows verified : {total_rows:,}")
        print(
            f"Columns       : "
            f"{len(header):,} -> "
            f"{len(header) - len(columns_to_remove):,}"
        )
        print(
            f"Removed       : "
            f"{', '.join(map(str, columns_to_remove))}"
        )
        print()
        print(
            "Verification result: only the requested columns were "
            "removed; row count, row order, remaining column order, "
            "and all remaining values are identical."
        )

    except KeyboardInterrupt:
        print("\n\nOperation cancelled by user.", file=sys.stderr)

        if temp_output_path.exists():
            temp_output_path.unlink()

        sys.exit(130)

    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)

        # Never leave an unverified temporary file behind.
        if temp_output_path.exists():
            try:
                temp_output_path.unlink()
            except OSError:
                pass

        sys.exit(1)


if __name__ == "__main__":
    main()