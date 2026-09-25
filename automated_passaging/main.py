# Copyright (c) 2026 Flavia Araujo
# SPDX-License-Identifier: MIT

import os
import sys
import watcher
import argparse
from pathlib import Path


def main():
    # Argument parsing from Momentum call
    parser = argparse.ArgumentParser(description="Watch a file and run a command when it changes.")
    parser.add_argument("target_OD", help="Target OD value for saturation (e.g., 3)")
    parser.add_argument("wells_threshold", help="Threshold percentage of wells (e.g., 70 for 70%)")
    parser.add_argument("target_volume", help="Target volume to dispense (e.g., 125)")
    parser.add_argument("process_arg", help="Momentum process to call when the threshold is met")
    args = parser.parse_args()

    # Request path interactively
    path_arg = input(" --- Enter path to the CSV file to watch (e.g., C://GP_20250930_111331_OD.csv): ").strip()
    if not path_arg:
        print("\n[GP Watcher] Error: no file path provided", file=sys.stderr)
        return 2

    target_path = Path(os.path.expanduser(path_arg))
    if target_path.is_dir():
        print("\n[GP Watcher] Error: path must be a file", file=sys.stderr)
        return 2

    parent_dir = target_path.resolve().parent
    if not parent_dir.exists():
        print(f"\n[GP Watcher] Error: directory does not exist: {parent_dir}", file=sys.stderr)
        return 2

    # target_path = ("/Users/flavia/PycharmProjects/growth_profile_watcher/sample/GP_20171118_012019_OD_24.csv")
    print(f"\n[GP Watcher] Watching: {target_path}")
    watcher.start_watching(
        target_path,
        args.process_arg,
        args.target_OD,
        args.wells_threshold,
        args.target_volume
    )


if __name__ == "__main__":
    raise SystemExit(main())

