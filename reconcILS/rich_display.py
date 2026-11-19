"""Rich terminal presentation helpers used by the reconcILS wrappers.

The reconciliation algorithm and file formats do not depend on Rich.  These
helpers provide DAFT-style panels, tables, and live progress when verbose mode
is enabled, with a plain-text fallback when Rich is unavailable.
"""
from __future__ import annotations

from contextlib import nullcontext
from pathlib import Path
import time
from typing import Any, Iterable, Sequence

try:
    from rich.console import Console, Group
    from rich.live import Live
    from rich.panel import Panel
    from rich.progress import (
        BarColumn,
        Progress,
        SpinnerColumn,
        TaskProgressColumn,
        TextColumn,
        TimeElapsedColumn,
        TimeRemainingColumn,
    )
    from rich.table import Table

    RICH_AVAILABLE = True
except Exception:  # pragma: no cover - exercised only in minimal environments
    Console = Group = Live = Panel = Progress = Table = None  # type: ignore[assignment]
    RICH_AVAILABLE = False


def format_seconds(seconds: float) -> str:
    """Format a duration as HH:MM:SS."""
    seconds = int(max(seconds, 0))
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    remaining = seconds % 60
    return f"{hours:02d}:{minutes:02d}:{remaining:02d}"


class RunDisplay:
    """Small DAFT-style panel/table facade with a plain-text fallback."""

    def __init__(self, enabled: bool = False, console: Any | None = None):
        self.enabled = bool(enabled)
        self.console = console
        if self.enabled and self.console is None and RICH_AVAILABLE:
            self.console = Console()

    @property
    def rich(self) -> bool:
        return self.enabled and self.console is not None and RICH_AVAILABLE

    def panel(self, title: str, message: str, style: str = "green") -> None:
        if not self.enabled:
            return
        if self.rich:
            self.console.print(Panel(message, title=title, border_style=style))
            return
        print()
        print("=" * 80)
        print(title)
        print("-" * 80)
        print(message)
        print("=" * 80)

    def table(
        self,
        title: str,
        columns: Sequence[str],
        rows: Iterable[Sequence[Any]],
        *,
        status_column: int | None = None,
    ) -> None:
        if not self.enabled:
            return
        materialized = [list(row) for row in rows]
        if self.rich:
            table = Table(title=title)
            for column in columns:
                table.add_column(str(column))
            for row in materialized:
                style = None
                if status_column is not None and status_column < len(row):
                    status = str(row[status_column]).upper()
                    if status == "OK":
                        style = "green"
                    elif status == "WARNING":
                        style = "yellow"
                    elif status in {"FAIL", "FAILED", "ERROR", "TIMEOUT"}:
                        style = "red"
                values = [str(value) for value in row]
                table.add_row(*values, style=style)
            self.console.print(table)
            return

        print()
        print(title)
        print("\t".join(str(column) for column in columns))
        for row in materialized:
            print("\t".join(str(value) for value in row))

    def message(self, message: str, style: str | None = None) -> None:
        if self.rich:
            self.console.print(message, style=style)
        else:
            print(message, flush=True)

    def status(self, message: str):
        if self.rich:
            return self.console.status(message, spinner="dots")
        return nullcontext()


class ParallelDashboard:
    """DAFT-style live progress dashboard for the parameterized batch example."""

    def __init__(
        self,
        *,
        total: int,
        workers: int,
        output: Path,
        state: Path,
        enabled: bool,
    ) -> None:
        self.total = total
        self.workers = workers
        self.output = output
        self.state = state
        self.display = RunDisplay(enabled)
        self.start_time = time.perf_counter()
        self.progress = None
        self.task_id = None
        self.live = None
        self.closed = False
        self.completed = 0
        self.saved = 0
        self.active = 0
        self.failed = 0
        self.checkpoint_rows = 0

        if self.display.rich:
            self.progress = Progress(
                SpinnerColumn(),
                TextColumn("[bold cyan]Reconciling gene trees"),
                BarColumn(),
                TaskProgressColumn(),
                TextColumn("[bold]{task.completed}/{task.total}"),
                TimeElapsedColumn(),
                TimeRemainingColumn(),
                expand=True,
            )
            self.task_id = self.progress.add_task("reconcILS", total=total)
            self.live = Live(
                self._render(),
                console=self.display.console,
                refresh_per_second=4,
                transient=False,
            )

    @property
    def rich(self) -> bool:
        return self.display.rich

    def _render(self):
        assert self.progress is not None
        assert Table is not None and Panel is not None and Group is not None
        table = Table.grid(padding=(0, 2))
        table.add_column(justify="left", style="bold")
        table.add_column(justify="right")
        table.add_row("Total gene trees", str(self.total))
        table.add_row("Processed", f"{self.completed}/{self.total}")
        table.add_row("Saved successfully", str(self.saved))
        table.add_row("Active workers", f"{self.active}/{self.workers}")
        table.add_row("Failures", str(self.failed))
        table.add_row("Checkpoint rows", str(self.checkpoint_rows))
        table.add_row("Elapsed wall time", format_seconds(time.perf_counter() - self.start_time))
        table.add_row("Output", str(self.output))
        table.add_row("State directory", str(self.state))
        return Panel.fit(
            Group(self.progress, table),
            title="[bold green]reconcILS parallel reconciliation[/bold green]",
            border_style="red",
            padding=(1, 2),
        )

    def start(self, *, completed: int = 0, saved: int = 0, active: int = 0) -> None:
        self.completed = completed
        self.saved = saved
        self.active = active
        self.checkpoint_rows = saved
        if self.rich:
            assert self.progress is not None and self.task_id is not None and self.live is not None
            self.progress.update(self.task_id, completed=completed)
            self.live.start()
            self.live.update(self._render())
        else:
            print(
                f"{self.total} selected trees; {saved} saved; {self.workers} workers; "
                f"output={self.output}; state={self.state}",
                flush=True,
            )

    def update(
        self,
        *,
        completed: int,
        saved: int,
        active: int,
        failed: int = 0,
        checkpoint_rows: int | None = None,
    ) -> None:
        self.completed = completed
        self.saved = saved
        self.active = active
        self.failed = failed
        if checkpoint_rows is not None:
            self.checkpoint_rows = checkpoint_rows
        if self.rich:
            assert self.progress is not None and self.task_id is not None and self.live is not None
            self.progress.update(self.task_id, completed=min(completed, self.total))
            self.live.update(self._render())

    def log(self, message: str, style: str | None = None) -> None:
        if self.rich:
            assert self.display.console is not None
            self.display.console.print(message, style=style)
        else:
            print(message, flush=True)

    def close(self) -> None:
        if self.closed:
            return
        self.closed = True
        if self.live is not None:
            self.live.stop()

    def finish(self, *, success: bool, message: str, warning: bool = False) -> None:
        self.close()
        if self.display.enabled:
            if self.display.rich:
                self.display.console.print()
            style = "yellow" if warning else ("green" if success else "red")
            title = (
                "reconcILS run interrupted"
                if warning
                else ("reconcILS run complete" if success else "reconcILS run failed")
            )
            self.display.panel(title, message, style)
        else:
            print(message, flush=True)
