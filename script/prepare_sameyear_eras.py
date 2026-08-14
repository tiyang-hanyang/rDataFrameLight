#!/usr/bin/env python3
import argparse
import json
import re
from pathlib import Path


YEAR_ERA_ORDER = {
    "2022": ["C", "D", "E", "F", "G"],
    "2023": ["C", "D"],
    "2024": ["C", "D", "E", "F", "G", "H", "I"],
}

ERA_PATTERN = re.compile(r"Run(20\d{2})([A-Z])")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Clone a single-era json into the other eras of the same year in the same directory. "
            "Only exact Run20XXY tokens are replaced."
        )
    )
    parser.add_argument(
        "--input-json",
        required=True,
        help="Template json path, typically a single-era file such as ...Run2024C.json",
    )
    parser.add_argument(
        "--target-eras",
        default=None,
        help=(
            "Optional comma-separated target eras, for example Run2024D,Run2024E. "
            "Default: all other eras in the same year after detecting the source era."
        ),
    )
    parser.add_argument(
        "--source-era",
        default=None,
        help="Optional explicit source era override, for example Run2024C.",
    )
    parser.add_argument(
        "--no-overwrite",
        action="store_true",
        help="Do not overwrite existing target files.",
    )
    return parser.parse_args()


def load_json(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def dump_json(path, payload):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=4)
        handle.write("\n")


def detect_source_era(input_path, payload, explicit_source_era=None):
    if explicit_source_era is not None:
        return explicit_source_era

    filename_match = ERA_PATTERN.search(input_path.name)
    if filename_match:
        return filename_match.group(0)

    era_field = payload.get("era")
    if isinstance(era_field, str):
        era_match = ERA_PATTERN.fullmatch(era_field)
        if era_match:
            return era_field
    if isinstance(era_field, list) and len(era_field) == 1 and isinstance(era_field[0], str):
        era_match = ERA_PATTERN.fullmatch(era_field[0])
        if era_match:
            return era_field[0]

    raise ValueError(
        "Failed to detect source era from filename or json content. Please pass --source-era explicitly."
    )


def build_default_targets(source_era):
    match = ERA_PATTERN.fullmatch(source_era)
    if not match:
        raise ValueError(f"Unsupported source era format: {source_era}")

    year = match.group(1)
    letter = match.group(2)
    if year not in YEAR_ERA_ORDER:
        raise ValueError(f"Unsupported year for same-year cloning: {year}")
    if letter not in YEAR_ERA_ORDER[year]:
        raise ValueError(f"Unsupported era letter {letter} for year {year}")

    return [f"Run{year}{target_letter}" for target_letter in YEAR_ERA_ORDER[year] if target_letter != letter]


def parse_target_eras(raw_target_eras, source_era):
    if raw_target_eras is None:
        return build_default_targets(source_era)

    target_eras = [item.strip() for item in raw_target_eras.split(",") if item.strip()]
    if not target_eras:
        raise ValueError("Empty --target-eras was provided.")

    source_year = ERA_PATTERN.fullmatch(source_era).group(1)
    for target_era in target_eras:
        target_match = ERA_PATTERN.fullmatch(target_era)
        if target_match is None:
            raise ValueError(f"Invalid target era format: {target_era}")
        if target_match.group(1) != source_year:
            raise ValueError(
                f"Target era {target_era} is not in the same year as source era {source_era}."
            )
        if target_era == source_era:
            raise ValueError(f"Target era list should not include the source era itself: {source_era}")

    return target_eras


def replace_era_token(value, source_era, target_era):
    if isinstance(value, str):
        return value.replace(source_era, target_era)
    if isinstance(value, list):
        return [replace_era_token(item, source_era, target_era) for item in value]
    if isinstance(value, dict):
        return {
            replace_era_token(key, source_era, target_era): replace_era_token(item, source_era, target_era)
            for key, item in value.items()
        }
    return value


def build_output_path(input_path, source_era, target_era):
    if source_era in input_path.name:
        return input_path.with_name(input_path.name.replace(source_era, target_era))
    return input_path.with_name(f"{input_path.stem}_{target_era}{input_path.suffix}")


def main():
    args = parse_args()
    input_path = Path(args.input_json).resolve()
    payload = load_json(input_path)

    source_era = detect_source_era(input_path, payload, args.source_era)
    target_eras = parse_target_eras(args.target_eras, source_era)
    overwrite = not args.no_overwrite

    print(f"input json: {input_path}")
    print(f"source era: {source_era}")
    print(f"overwrite: {overwrite}")

    for target_era in target_eras:
        output_path = build_output_path(input_path, source_era, target_era)
        existed_before = output_path.exists()
        if existed_before and not overwrite:
            print(f"skip existing: {output_path}")
            continue

        output_payload = replace_era_token(payload, source_era, target_era)
        dump_json(output_path, output_payload)
        action = "overwrite" if existed_before else "write"
        print(f"{action}: {output_path}")


if __name__ == "__main__":
    main()
