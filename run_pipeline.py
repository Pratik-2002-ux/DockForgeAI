#!/usr/bin/env python3
"""
DockForgeAI one-click pipeline runner.

Usage:
    python run_pipeline.py --config config/pipeline.yaml
    python run_pipeline.py --config config/pipeline.yaml --dry-run
    python run_pipeline.py --config config/pipeline.yaml --resume
    python run_pipeline.py --config config/pipeline.yaml --from-stage docking
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:
    sys.exit(
        "PyYAML is not installed.\n"
        "Install it with:\n"
        "    python -m pip install pyyaml"
    )


def load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {path}")

    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)

    if not isinstance(data, dict):
        raise ValueError("The YAML configuration must contain a mapping.")

    return data


def safe_name(value: str) -> str:
    return "".join(
        character if character.isalnum() or character in "-_" else "_"
        for character in value
    )


def command_program(command: str) -> str | None:
    try:
        parts = shlex.split(command)
    except ValueError:
        return None

    return parts[0] if parts else None


def command_exists(command: str, project_root: Path) -> bool:
    program = command_program(command)

    if not program:
        return False

    if "/" in program:
        candidate = Path(program)
        if not candidate.is_absolute():
            candidate = project_root / candidate
        return candidate.exists()

    return shutil.which(program) is not None


def expand_command(
    command: str,
    variables: dict[str, Any],
) -> str:
    expanded = command

    for key, value in variables.items():
        expanded = expanded.replace(
            "{" + str(key) + "}",
            str(value),
        )

    return os.path.expandvars(os.path.expanduser(expanded))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def run_command(
    command: str,
    working_directory: Path,
    log_path: Path,
    environment: dict[str, str],
    dry_run: bool,
) -> int:
    print(f"\n$ {command}")

    if dry_run:
        return 0

    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("w", encoding="utf-8") as log_handle:
        process = subprocess.Popen(
            command,
            cwd=working_directory,
            shell=True,
            executable="/bin/zsh",
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=environment,
        )

        assert process.stdout is not None

        for line in process.stdout:
            print(line, end="")
            log_handle.write(line)
            log_handle.flush()

        return process.wait()


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run DockForgeAI stages from a YAML configuration."
    )

    parser.add_argument(
        "--config",
        required=True,
        help="Path to pipeline YAML configuration.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Display commands without executing them.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip stages already marked as completed.",
    )
    parser.add_argument(
        "--from-stage",
        default=None,
        help="Start from this stage name.",
    )
    parser.add_argument(
        "--to-stage",
        default=None,
        help="Stop after this stage name.",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_arguments()

    config_path = Path(args.config).resolve()
    config = load_yaml(config_path)

    project_root = Path(
        config.get("project_root", config_path.parent.parent)
    ).expanduser().resolve()

    if not project_root.exists():
        raise FileNotFoundError(
            f"Project root does not exist: {project_root}"
        )

    pipeline_name = str(
        config.get("pipeline_name", "DockForgeAI")
    )

    variables = dict(config.get("variables", {}))
    variables["project_root"] = project_root
    variables["config_path"] = config_path

    output_root = Path(
        expand_command(
            str(
                config.get(
                    "output_root",
                    project_root / "runs",
                )
            ),
            variables,
        )
    )

    if not output_root.is_absolute():
        output_root = project_root / output_root

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_directory = output_root / run_id
    logs_directory = run_directory / "logs"
    state_path = output_root / "pipeline_state.json"

    output_root.mkdir(parents=True, exist_ok=True)
    run_directory.mkdir(parents=True, exist_ok=True)
    logs_directory.mkdir(parents=True, exist_ok=True)

    environment = os.environ.copy()

    for key, value in config.get("environment", {}).items():
        environment[str(key)] = expand_command(
            str(value),
            variables,
        )

    stages = config.get("stages")

    if not isinstance(stages, list) or not stages:
        raise ValueError(
            "The configuration must contain a non-empty 'stages' list."
        )

    previous_state: dict[str, Any] = {}

    if state_path.exists():
        try:
            previous_state = json.loads(
                state_path.read_text(encoding="utf-8")
            )
        except json.JSONDecodeError:
            previous_state = {}

    completed_stages = set(
        previous_state.get("completed_stages", [])
    )

    selected_stages: list[dict[str, Any]] = []
    started = args.from_stage is None

    for stage in stages:
        if not isinstance(stage, dict):
            raise ValueError("Every stage must be a YAML mapping.")

        stage_name = str(stage.get("name", "")).strip()

        if not stage_name:
            raise ValueError("Each stage requires a name.")

        if args.from_stage and stage_name == args.from_stage:
            started = True

        if not started:
            continue

        selected_stages.append(stage)

        if args.to_stage and stage_name == args.to_stage:
            break

    if args.from_stage and not selected_stages:
        raise ValueError(
            f"Starting stage was not found: {args.from_stage}"
        )

    manifest: dict[str, Any] = {
        "pipeline_name": pipeline_name,
        "run_id": run_id,
        "started_at": datetime.now().isoformat(),
        "project_root": str(project_root),
        "config": str(config_path),
        "dry_run": args.dry_run,
        "resume": args.resume,
        "stages": [],
    }

    print("=" * 72)
    print(f"{pipeline_name}")
    print(f"Project root : {project_root}")
    print(f"Run ID       : {run_id}")
    print(f"Logs         : {logs_directory}")
    print("=" * 72)

    for position, stage in enumerate(selected_stages, start=1):
        stage_name = str(stage["name"])
        enabled = bool(stage.get("enabled", True))

        if not enabled:
            print(f"\n[SKIP] {stage_name}: disabled")
            continue

        if args.resume and stage_name in completed_stages:
            print(f"\n[SKIP] {stage_name}: already completed")
            continue

        commands = stage.get("commands", [])

        if isinstance(commands, str):
            commands = [commands]

        if not isinstance(commands, list) or not commands:
            raise ValueError(
                f"Stage '{stage_name}' has no commands."
            )

        working_directory_value = expand_command(
            str(stage.get("working_directory", project_root)),
            variables,
        )
        working_directory = Path(working_directory_value)

        if not working_directory.is_absolute():
            working_directory = project_root / working_directory

        if not working_directory.exists():
            raise FileNotFoundError(
                f"Working directory for '{stage_name}' does not exist: "
                f"{working_directory}"
            )

        print(
            f"\n[{position}/{len(selected_stages)}] "
            f"Starting: {stage_name}"
        )

        stage_started = time.time()
        stage_record: dict[str, Any] = {
            "name": stage_name,
            "started_at": datetime.now().isoformat(),
            "commands": [],
        }

        for command_number, raw_command in enumerate(
            commands,
            start=1,
        ):
            command = expand_command(
                str(raw_command),
                variables,
            )

            if not command_exists(command, project_root):
                program = command_program(command)
                raise FileNotFoundError(
                    f"Required command is unavailable for stage "
                    f"'{stage_name}': {program}"
                )

            log_name = (
                f"{position:02d}_"
                f"{safe_name(stage_name)}_"
                f"{command_number:02d}.log"
            )
            log_path = logs_directory / log_name

            command_started = time.time()
            return_code = run_command(
                command=command,
                working_directory=working_directory,
                log_path=log_path,
                environment=environment,
                dry_run=args.dry_run,
            )

            command_record = {
                "command": command,
                "return_code": return_code,
                "log": str(log_path),
                "duration_seconds": round(
                    time.time() - command_started,
                    3,
                ),
            }
            stage_record["commands"].append(command_record)

            if return_code != 0:
                stage_record["status"] = "failed"
                stage_record["finished_at"] = (
                    datetime.now().isoformat()
                )
                manifest["stages"].append(stage_record)
                manifest["status"] = "failed"
                manifest["finished_at"] = (
                    datetime.now().isoformat()
                )

                write_json(
                    run_directory / "run_manifest.json",
                    manifest,
                )

                print(
                    f"\n[FAILED] Stage '{stage_name}' "
                    f"returned exit code {return_code}."
                )
                print(f"Log: {log_path}")
                return return_code

        stage_record["status"] = (
            "dry_run" if args.dry_run else "complete"
        )
        stage_record["finished_at"] = datetime.now().isoformat()
        stage_record["duration_seconds"] = round(
            time.time() - stage_started,
            3,
        )
        manifest["stages"].append(stage_record)

        if not args.dry_run:
            completed_stages.add(stage_name)

            write_json(
                state_path,
                {
                    "pipeline_name": pipeline_name,
                    "completed_stages": sorted(completed_stages),
                    "updated_at": datetime.now().isoformat(),
                },
            )

        print(f"[DONE] {stage_name}")

    manifest["status"] = (
        "dry_run" if args.dry_run else "complete"
    )
    manifest["finished_at"] = datetime.now().isoformat()

    write_json(
        run_directory / "run_manifest.json",
        manifest,
    )

    print("\n" + "=" * 72)
    print("DockForgeAI pipeline finished successfully.")
    print(f"Run manifest: {run_directory / 'run_manifest.json'}")
    print("=" * 72)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nPipeline interrupted by the user.")
        raise SystemExit(130)
    except Exception as error:
        print(f"\nPipeline error: {error}", file=sys.stderr)
        raise SystemExit(1)
