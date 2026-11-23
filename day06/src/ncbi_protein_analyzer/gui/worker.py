"""Thread-safe background task runner bridged to Tkinter via queue + root.after."""

from __future__ import annotations

import queue
import threading
import tkinter as tk
from collections.abc import Callable
from typing import Any

from ncbi_protein_analyzer.ncbi.client import CancelledError

ProgressFn = Callable[[str, float | None], None]
TaskFn = Callable[[threading.Event, ProgressFn], Any]
SuccessCb = Callable[[Any], None]
ErrorCb = Callable[[Exception], None]
CancelCb = Callable[[str], None]


class BackgroundRunner:
    """Runs blocking tasks on a daemon thread and delivers updates on the Tk main thread."""

    def __init__(self, root: tk.Misc, poll_interval_ms: int = 80) -> None:
        self._root = root
        self._poll_interval_ms = poll_interval_ms
        self._queue: queue.Queue[tuple[str, Any]] = queue.Queue()
        self._cancel_event: threading.Event | None = None
        self._thread: threading.Thread | None = None
        self._on_progress: ProgressFn | None = None
        self._on_success: SuccessCb | None = None
        self._on_error: ErrorCb | None = None
        self._on_cancel: CancelCb | None = None

    @property
    def is_running(self) -> bool:
        return self._thread is not None and self._thread.is_alive()

    def start(
        self,
        task: TaskFn,
        *,
        on_progress: ProgressFn,
        on_success: SuccessCb,
        on_error: ErrorCb,
        on_cancel: CancelCb,
    ) -> bool:
        """Start a background task if none is currently running."""
        if self.is_running:
            return False

        self._queue = queue.Queue()
        self._cancel_event = threading.Event()
        self._on_progress = on_progress
        self._on_success = on_success
        self._on_error = on_error
        self._on_cancel = on_cancel

        cancel_event = self._cancel_event

        def _emit_progress(message: str, fraction: float | None = None) -> None:
            self._queue.put(("progress", (message, fraction)))

        def _worker() -> None:
            try:
                result = task(cancel_event, _emit_progress)
            except CancelledError as exc:
                self._queue.put(("cancelled", str(exc)))
            except Exception as exc:  # noqa: BLE001 - surface any failure cleanly to GUI
                self._queue.put(("error", exc))
            else:
                self._queue.put(("done", result))

        self._thread = threading.Thread(target=_worker, daemon=True)
        self._thread.start()
        self._root.after(self._poll_interval_ms, self._poll)
        return True

    def cancel(self) -> None:
        """Signal the active background task to cancel."""
        if self._cancel_event is not None:
            self._cancel_event.set()

    def _poll(self) -> None:
        finished = False
        while True:
            try:
                kind, payload = self._queue.get_nowait()
            except queue.Empty:
                break

            if kind == "progress" and self._on_progress is not None:
                message, fraction = payload
                self._on_progress(message, fraction)
            elif kind == "done":
                finished = True
                if self._on_success is not None:
                    self._on_success(payload)
            elif kind == "cancelled":
                finished = True
                if self._on_cancel is not None:
                    self._on_cancel(str(payload))
            elif kind == "error":
                finished = True
                if self._on_error is not None:
                    self._on_error(payload)

        if not finished and self.is_running:
            self._root.after(self._poll_interval_ms, self._poll)
        elif finished:
            self._thread = None
            self._cancel_event = None
