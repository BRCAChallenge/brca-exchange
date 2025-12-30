#!/usr/bin/env python3
"""
Generate a new BRCA Exchange data release.

This script sets up the environment and kicks off the pipeline to create
a new data release. It handles:
- Creating a working directory
- Cloning the code repository
- Generating the pipeline configuration from template
- Running the build-release Make target

Python 3.13+ compatible.
"""

import argparse
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

try:
    from jinja2 import Environment, FileSystemLoader, select_autoescape
except ImportError:
    print("Error: jinja2 is required. Install with: pip install jinja2", file=sys.stderr)
    sys.exit(1)


def resolve_path(path: str) -> Path:
    """Resolve a path to its absolute form."""
    return Path(path).resolve()


def run_command(
    cmd: list[str],
    cwd: Optional[Path] = None,
    check: bool = True,
    capture_output: bool = False
) -> subprocess.CompletedProcess:
    """
    Run a shell command with proper error handling.

    Args:
        cmd: Command and arguments as a list
        cwd: Working directory for the command
        check: Whether to raise exception on non-zero exit
        capture_output: Whether to capture stdout/stderr

    Returns:
        CompletedProcess instance
    """
    print(f"Running: {' '.join(cmd)}")
    return subprocess.run(
        cmd,
        cwd=cwd,
        check=check,
        capture_output=capture_output,
        text=True
    )


def clone_or_update_repo(code_base: Path, git_commit: str) -> None:
    """
    Clone the repository if it doesn't exist, then checkout the specified commit.

    Args:
        code_base: Path where the code should be cloned
        git_commit: Git commit/branch/tag to checkout
    """
    repo_url = "https://github.com/BRCAChallenge/brca-exchange-kb.git"

    if not code_base.exists():
        print(f"Cloning repository to {code_base}...")
        run_command(["git", "clone", repo_url, str(code_base)])
    else:
        print(f"Repository already exists at {code_base}")
        print("Fetching latest changes from remote...")
        run_command(["git", "fetch", "origin"], cwd=code_base)

    print(f"Checking out {git_commit}...")
    # Try to checkout directly first
    try:
        run_command(["git", "checkout", git_commit], cwd=code_base)
    except subprocess.CalledProcessError:
        # If direct checkout fails, try as a remote branch
        print(f"Direct checkout failed, trying origin/{git_commit}...")
        run_command(["git", "checkout", "-b", git_commit, f"origin/{git_commit}"], cwd=code_base)


def generate_config(
    template_path: Path,
    output_path: Path,
    context: dict[str, str]
) -> None:
    """
    Generate configuration file from Jinja2 template.

    Args:
        template_path: Path to the Jinja2 template file
        output_path: Path where the generated config should be written
        context: Dictionary of template variables
    """
    print(f"Generating configuration file: {output_path}")

    # Set up Jinja2 environment
    template_dir = template_path.parent
    env = Environment(
        loader=FileSystemLoader(str(template_dir)),
        autoescape=select_autoescape(),
        keep_trailing_newline=True
    )

    # Load and render template
    template = env.get_template(template_path.name)
    rendered = template.render(**context)

    # Write output
    output_path.write_text(rendered)
    print(f"Configuration written to {output_path}")


def run_pipeline(
    code_base: Path,
    gene_config_filename: str
) -> None:
    """
    Execute the pipeline build-release target.

    Args:
        code_base: Path to the code repository
        gene_config_filename: Name of the gene configuration file
    """
    pipeline_dir = code_base / "pipeline"

    print("\nKicking off pipeline!")
    print(f"Gene configuration: {gene_config_filename}")

    run_command(
        ["make", f"GENE_CONFIG_FILENAME={gene_config_filename}", "build-release"],
        cwd=pipeline_dir
    )


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Generate a new BRCA Exchange data release",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  %(prog)s /data/releases /data/credentials /data/previous_releases
  %(prog)s /data/releases /data/credentials /data/previous_releases gene_config_brca_hbop.txt v1.2.3
  %(prog)s /data/releases /data/credentials /data/previous_releases gene_config_brca_only.txt my-feature-branch
        """
    )

    parser.add_argument(
        "root_dir",
        type=str,
        help="Root directory for the release (working directory will be created here)"
    )

    parser.add_argument(
        "credentials_path",
        type=str,
        help="Path to Luigi credentials configuration file"
    )

    parser.add_argument(
        "previous_release_dir",
        type=str,
        help="Directory containing previous release for comparison"
    )

    parser.add_argument(
        "gene_config_filename",
        type=str,
        nargs='?',
        default="gene_config_brca_only.txt",
        help="Gene configuration filename (default: gene_config_brca_only.txt)"
    )

    parser.add_argument(
        "git_commit",
        type=str,
        nargs='?',
        default="master",
        help="Git commit/branch/tag to checkout from GitHub (default: master)"
    )

    args = parser.parse_args()

    # Resolve all paths to absolute
    root_dir = resolve_path(args.root_dir)
    credentials_path = resolve_path(args.credentials_path)
    previous_release_dir = resolve_path(args.previous_release_dir)

    # Generate data date and working directory
    data_date = datetime.now().strftime("%Y-%m-%d")
    work_dir = root_dir / f"data_release_{data_date}"

    print(f"=== BRCA Exchange Release Generator ===")
    print(f"Data Date: {data_date}")
    print(f"Working Directory: {work_dir}")
    print(f"Gene Configuration: {args.gene_config_filename}")
    print(f"Git Commit: {args.git_commit}")
    print("=" * 40)

    # Create working directory
    work_dir.mkdir(parents=True, exist_ok=True)
    print(f"Created working directory: {work_dir}")

    # Set up code base
    code_base = work_dir / "code"
    clone_or_update_repo(code_base, args.git_commit)

    # Prepare template context
    template_path = code_base / "pipeline" / "pipeline_running" / "brca_pipeline_cfg.mk.j2"
    config_path = code_base / "pipeline" / "brca_pipeline_cfg.mk"

    context = {
        "DATA_DATE": data_date,
        "WORK_DIR": str(work_dir),
        "CODE_BASE": str(code_base),
        "CREDENTIALS_PATH": str(credentials_path),
        "PREVIOUS_RELEASE_DIR": str(previous_release_dir),
        "GIT_COMMIT": args.git_commit,
    }

    # Generate configuration
    generate_config(template_path, config_path, context)

    # Print usage information
    print("\n" + "=" * 40)
    print("Configuration generated successfully!")
    print("\nYou can issue pipeline commands using:")
    print(f"  make CONFIG_PATH={config_path} [cmd]")
    print("-- or --")
    print(f"  cd {code_base / 'pipeline'} && make [cmd]")
    print("=" * 40)

    # Run the pipeline
    try:
        run_pipeline(code_base, args.gene_config_filename)
        print("\n" + "=" * 40)
        print("Pipeline completed successfully!")
        print("=" * 40)
        return 0
    except subprocess.CalledProcessError as e:
        print(f"\nError: Pipeline failed with exit code {e.returncode}", file=sys.stderr)
        return e.returncode
    except KeyboardInterrupt:
        print("\n\nInterrupted by user", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
