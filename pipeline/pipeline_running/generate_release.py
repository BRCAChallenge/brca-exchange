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

try:
    import psycopg2
    from psycopg2 import sql
except ImportError:
    print("Error: psycopg2 is required. Install with: pip install psycopg2", file=sys.stderr)
    sys.exit(1)


# Same default/override convention as the rest of the pipeline's DB-touching
# scripts (e.g. variant_analysis/run_vep_analysis.py) -- a single Postgres
# instance holding every release's tables, disambiguated by schema.
DEFAULT_PIPELINE_DB_URL = "postgresql://postgres:postgres@localhost/storage.pg"


def resolve_path(path: str) -> Path:
    """Resolve a path to its absolute form."""
    return Path(path).resolve()


def get_current_branch(script_dir: Path) -> str:
    """
    Determine the git branch of the checkout this script is running from,
    so running from a feature branch deploys that branch by default rather
    than silently falling back to master.

    Falls back to "master" if the branch can't be determined (e.g. the
    script isn't inside a git checkout, or HEAD is detached).
    """
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            cwd=script_dir,
            check=True,
            capture_output=True,
            text=True,
        )
        branch = result.stdout.strip()
        if branch and branch != "HEAD":
            return branch
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass
    return "master"


def run_command(
    cmd: list[str],
    cwd: Optional[Path] = None,
    check: bool = True,
    capture_output: bool = False,
    env: Optional[dict[str, str]] = None,
) -> subprocess.CompletedProcess:
    """
    Run a shell command with proper error handling.

    Args:
        cmd: Command and arguments as a list
        cwd: Working directory for the command
        check: Whether to raise exception on non-zero exit
        capture_output: Whether to capture stdout/stderr
        env: Environment for the subprocess (default: inherit this
            process's environment)

    Returns:
        CompletedProcess instance
    """
    print(f"Running: {' '.join(cmd)}")
    return subprocess.run(
        cmd,
        cwd=cwd,
        check=check,
        capture_output=capture_output,
        text=True,
        env=env,
    )


def clone_or_update_repo(code_base: Path, git_commit: str) -> None:
    """
    Clone the repository if it doesn't exist, then check out the specified
    branch/commit/tag and fast-forward it to match origin.

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
        return

    # `git checkout` above is a no-op if code_base was already on this
    # branch, so a reused work directory would otherwise keep running
    # whatever commit it happened to be on -- even after the `git fetch`
    # above -- silently ignoring newer commits on origin. Fast-forward
    # explicitly to pick those up.
    try:
        run_command(["git", "merge", "--ff-only", f"origin/{git_commit}"], cwd=code_base)
    except subprocess.CalledProcessError:
        print(f"Note: no origin/{git_commit} to fast-forward from "
              "(may be a tag or a specific commit) -- leaving as checked out.")


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


def create_schema_namespace(db_url: str, schema_name: str) -> None:
    """
    Create this release's (empty) PostgreSQL schema namespace. Table
    structure is added later, by apply_migrations() -- Django's `data` app
    models are the source of truth for that DDL.

    Idempotent: CREATE SCHEMA IF NOT EXISTS, safe to rerun against a work
    directory whose schema already exists (e.g. a reused work directory).

    Args:
        db_url: Postgres connection string.
        schema_name: Name of the schema to create for this release.
    """
    print(f"Creating PostgreSQL schema '{schema_name}'...")

    conn = psycopg2.connect(db_url)
    try:
        conn.autocommit = True
        with conn.cursor() as cur:
            cur.execute(
                sql.SQL("CREATE SCHEMA IF NOT EXISTS {}")
                .format(sql.Identifier(schema_name))
            )
    finally:
        conn.close()

    print(f"Schema '{schema_name}' created.")


def apply_migrations(code_base: Path, schema_name: str) -> None:
    """
    Run the `data` app's Django migrations against this release's schema,
    via the 'pipeline' DB alias (django/brca/site_settings.py points its
    search_path at $PIPELINE_SCHEMA). This is what actually creates the
    tables -- on a brand-new schema this issues real DDL matching current
    models.py exactly; naturally idempotent if rerun against a work
    directory that's already been migrated.

    Args:
        code_base: Path to the code repository (must already be checked
            out -- this shells out to its django/manage.py).
        schema_name: Name of the schema to migrate.
    """
    print(f"Running migrations against schema '{schema_name}'...")
    django_dir = code_base / "django"
    run_command(
        [sys.executable, "manage.py", "migrate", "--database=pipeline", "data"],
        cwd=django_dir,
        env={**os.environ, "PIPELINE_SCHEMA": schema_name},
    )


def spawn_pipeline(code_base: Path, work_dir: Path) -> Path:
    """
    Launch the pipeline build-release target in the background and return
    immediately -- this does not wait for the pipeline to finish.

    The gene configuration is baked into the generated brca_pipeline_cfg.mk
    (via GENE_CONFIG_FILENAME in the template context), not passed as a make
    variable here -- the Makefile no longer reads a GENE_CONFIG_FILENAME
    override directly.

    Both logs are now written by the Makefile itself: run-pipeline's Luigi
    output goes to PIPELINE_LOG, and the output of build-releases's subcommands
    (checkout, resource downloads, docker service startup, ...) goes to
    BUILD_RELEASE_LOG. (The output of make build-release itself is lost to /dev/null.)

    Args:
        code_base: Path to the code repository
        work_dir: This release's working directory. PIPELINE_LOG lives
            directly under it; BUILD_RELEASE_LOG lives under its tmp/
            subdirectory.

    Returns:
        Path to the log file the running Luigi pipeline is writing to.
    """
    pipeline_dir = code_base / "pipeline"

    log_dir = work_dir / "tmp"
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    pipeline_log = work_dir / f"pipeline_{timestamp}.log"
    build_release_log = log_dir / f"build_release_{timestamp}.log"

    print(f"\nKicking off pipeline! (spawned in the background, not waiting)")
    cmd = [
        "make",
        f"PIPELINE_LOG={pipeline_log}",
        f"BUILD_RELEASE_LOG={build_release_log}",
        "build-release",
    ]
    print(f"Running: {' '.join(cmd)}")
    print(f"build-release log: {build_release_log}")
    # Both logs above are opened by the Makefile's own recipes, not here --
    # nothing meaningful is left on this process's stdout/stderr to capture.
    subprocess.Popen(
        cmd,
        cwd=pipeline_dir,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        stdin=subprocess.DEVNULL,
        start_new_session=True,
    )

    return pipeline_log


def build_arg_parser(git_commit: str) -> argparse.ArgumentParser:
    """Build the CLI parser. `git_commit` is only used for the help text."""
    parser = argparse.ArgumentParser(
        description="Generate a new BRCA Exchange data release",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Builds from the current branch of this checkout ('{git_commit}').

Example usage:
  %(prog)s /data/releases /data/credentials
  %(prog)s /data/releases /data/credentials --previous-release-dir /data/previous_releases
  %(prog)s /data/releases /data/credentials gene_config_brca_hbop.txt --previous-release-dir /data/previous_releases
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
        "--previous-release-dir",
        dest="previous_release_dir",
        type=str,
        default=None,
        help="Directory containing previous release for comparison "
             "(optional; defaults to the working directory itself, i.e. no separate comparison dir)"
    )

    parser.add_argument(
        "--db-url",
        dest="db_url",
        type=str,
        default=os.environ.get("PIPELINE_DB_URL", DEFAULT_PIPELINE_DB_URL),
        help="Postgres connection string the release's PIPELINE_DB_SCHEMA is created in "
             "(default: $PIPELINE_DB_URL, falling back to %(default)r)"
    )

    parser.add_argument(
        "--run-gnomad-coverage-estimation",
        dest="run_gnomad_coverage_estimation",
        action="store_true",
        help="Recompute the gnomAD coverage parquets from raw gnomAD coverage data "
             "(DownloadGnomADCoverage in workflow/variant_assembly.py) instead of the "
             "default of downloading the pre-built parquets from brcaexchange.org. "
             "Slow -- downloads and reprocesses full genome/exome coverage summaries; "
             "the resulting parquets are copied back into the resources directory for "
             "future releases to reuse. Default: off (download the pre-built parquets)."
    )

    parser.add_argument(
        "gene_config_filename",
        type=str,
        nargs='?',
        default="gene_config_brca_only.txt",
        help="Gene configuration filename (default: gene_config_brca_only.txt)"
    )

    return parser


def print_run_summary(
    data_date: str,
    work_dir: Path,
    gene_config_filename: str,
    git_commit: str,
    previous_release_dir: Optional[Path],
    run_gnomad_coverage_estimation: bool = False,
) -> None:
    print(f"=== BRCA Exchange Release Generator ===")
    print(f"Data Date: {data_date}")
    print(f"Working Directory: {work_dir}")
    print(f"PIPELINE_DB_SCHEMA: {work_dir.name}")
    print(f"Gene Configuration: {gene_config_filename}")
    print(f"Git Commit: {git_commit}")
    print(f"Previous Release Dir: {previous_release_dir or '(none -- defaulting to working directory)'}")
    print(f"Recompute gnomAD Coverage (vs. download pre-built): {run_gnomad_coverage_estimation}")
    print("=" * 40)


def build_template_context(
    data_date: str,
    work_dir: Path,
    code_base: Path,
    credentials_path: Path,
    git_commit: str,
    gene_config_filename: str,
    db_schema: str,
    previous_release_dir: Optional[Path],
    run_gnomad_coverage_estimation: bool = False,
) -> dict[str, str]:
    """Build the Jinja2 context for rendering brca_pipeline_cfg.mk.j2."""
    context = {
        "DATA_DATE": data_date,
        "WORK_DIR": str(work_dir),
        "CODE_BASE": str(code_base),
        "CREDENTIALS_PATH": str(credentials_path),
        "GIT_COMMIT": git_commit,
        "GENE_CONFIG_FILENAME": gene_config_filename,
        "DB_SCHEMA": db_schema,
        "RUN_GNOMAD_COVERAGE_ESTIMATION": "true" if run_gnomad_coverage_estimation else "false",
    }
    if previous_release_dir is not None:
        context["PREVIOUS_RELEASE_DIR"] = str(previous_release_dir)
    return context


def print_config_ready(config_path: Path, pipeline_dir: Path) -> None:
    print("\n" + "=" * 40)
    print("Configuration generated successfully!")
    print("\nYou can issue pipeline commands using:")
    print(f"  make CONFIG_PATH={config_path} [cmd]")
    print("-- or --")
    print(f"  cd {pipeline_dir} && make [cmd]")
    print("=" * 40)


def print_pipeline_launched(pipeline_log: Path) -> None:
    print("\n" + "=" * 40)
    print("Pipeline launched in the background (not waiting for it to finish).")
    print(f"Follow its progress with:\n  tail -f {pipeline_log}")
    print("=" * 40)


def main() -> int:
    """Main entry point."""
    # The branch to build always comes from the checkout this script is
    # itself running from -- no CLI override, so there's no way to
    # accidentally deploy a different branch than the one you're standing in.
    git_commit = get_current_branch(Path(__file__).resolve().parent)
    args = build_arg_parser(git_commit).parse_args()

    # Resolve all paths to absolute
    root_dir = resolve_path(args.root_dir)
    credentials_path = resolve_path(args.credentials_path)
    previous_release_dir = resolve_path(args.previous_release_dir) if args.previous_release_dir else None

    # Generate data date and working directory. The working directory name
    # doubles as the pipeline's PostgreSQL schema name (see
    # brca_pipeline_cfg.mk.j2's PIPELINE_DB_SCHEMA), so it's built from an
    # underscore-separated date rather than DATA_DATE's hyphenated ISO form --
    # unquoted Postgres identifiers can't contain hyphens.
    data_date = datetime.now().strftime("%Y-%m-%d")
    work_dir = root_dir / f"data_release_{data_date.replace('-', '_')}"
    db_schema = work_dir.name

    print_run_summary(data_date, work_dir, args.gene_config_filename, git_commit, previous_release_dir,
                       args.run_gnomad_coverage_estimation)

    # Create working directory
    work_dir.mkdir(parents=True, exist_ok=True)
    print(f"Created working directory: {work_dir}")

    # Create this release's (empty) PostgreSQL schema namespace (named
    # after the work directory, same as brca_pipeline_cfg.mk.j2's own
    # DB_SCHEMA default) before doing anything else, so a bad --db-url
    # fails fast.
    try:
        create_schema_namespace(args.db_url, db_schema)
    except psycopg2.Error as e:
        print(f"\nError: could not create schema '{db_schema}': {e}", file=sys.stderr)
        return 1

    # Set up code base
    code_base = work_dir / "code"
    clone_or_update_repo(code_base, git_commit)

    # Populate the schema via Django's `data` app migrations, now that
    # code_base/django/manage.py exists to run them with.
    apply_migrations(code_base, db_schema)

    # Generate configuration
    template_path = code_base / "pipeline" / "pipeline_running" / "brca_pipeline_cfg.mk.j2"
    config_path = code_base / "pipeline" / "brca_pipeline_cfg.mk"
    context = build_template_context(
        data_date, work_dir, code_base, credentials_path, git_commit,
        args.gene_config_filename, db_schema, previous_release_dir,
        args.run_gnomad_coverage_estimation,
    )
    generate_config(template_path, config_path, context)
    print_config_ready(config_path, code_base / "pipeline")

    # Launch the pipeline in the background. build-release (particularly the
    # actual Luigi run within it) can take hours, so we don't wait for it --
    # just report where to watch it.
    try:
        pipeline_log = spawn_pipeline(code_base, work_dir)
    except OSError as e:
        print(f"\nError: could not launch the pipeline: {e}", file=sys.stderr)
        return 1

    print_pipeline_launched(pipeline_log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
