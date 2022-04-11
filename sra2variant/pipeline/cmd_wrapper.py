import os
import shutil
import subprocess
import logging
from abc import ABC, abstractmethod

from sra2variant.artifacts.base_file import _FileArtifacts


logging.basicConfig(
    format="[%(asctime)s]: %(message)s",
    datefmt="%Y-%m-%d %I:%M:%S %p",
    level=logging.INFO
)


class CMDwrapperBase(ABC):

    exec_name: str = None
    threads: str = "2"

    __slots__ = ["input_files", "output_files",
                 "cmd", "stdout", "stderr"]

    def __init__(
        self,
        input_files: _FileArtifacts,
        *args: str
    ) -> None:
        self.input_files: _FileArtifacts = input_files
        self.cmd: tuple[str] = (
            *self.exec_name.split(" "),
            *tuple(
                fn if os.path.isabs(fn) or not os.path.exists(fn)
                else os.path.relpath(fn, input_files.cwd)
                for fn in args
            )
        )
        self.stdout: str = "No stdout"
        self.stderr: str = "No stderr"

    def __str__(self) -> str:
        return " ".join(self.cmd)

    def __repr__(self) -> str:
        return str(self)

    @classmethod
    def set_exec_cmd(cls, exec_name: str) -> None:
        cls.exec_name = exec_name

    @classmethod
    def set_threads(cls, threads: int) -> None:
        cls.threads = str(threads)

    def execute_cmd(self) -> _FileArtifacts:
        logging.info(f"{self.input_files.workding_id} {self.exec_name}")
        with subprocess.Popen(
            self.cmd,
            text=True,
            cwd=self.input_files.cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ) as p:
            self.stdout, self.stderr = p.communicate()
        with open(self.input_files.log_file(), "a") as f:
            f.write(f"[Command]:\n{str(self)}\n")
            f.write(f"[stdout]:\n{self.stdout}\n")
            f.write(f"[stderr]:\n{self.stderr}\n")
        if p.returncode == 0:
            self._post_execution()
        return self.output_files

    @abstractmethod
    def _post_execution(self) -> None:
        return NotImplemented


class ErrorTolerance:

    max_errors: int = 0

    __slots__ = ["error_dir", "task_log_file"]

    def __init__(self, error_dir: str, task_log_file: str) -> None:
        self.error_dir: str = error_dir
        self.task_log_file: str = task_log_file

    @classmethod
    def set_max_errors(cls, max_errors: int) -> None:
        cls.max_errors = max_errors

    def handle(self, e: Exception) -> None:
        with open(self.task_log_file, "a") as f:
            f.write(f"Raised error: {e}")
        copied_log_file = os.path.join(
            self.error_dir,
            os.path.basename(self.task_log_file)
        )
        shutil.copyfile(self.task_log_file, copied_log_file)
        n_errors = len(os.listdir(self.error_dir))
        if n_errors > self.max_errors:
            raise RuntimeError(f"{n_errors} errors occurred exceeding maximum")
