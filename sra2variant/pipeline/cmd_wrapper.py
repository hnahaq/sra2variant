import os
import shutil
import subprocess
import csv
from abc import ABC, abstractmethod


CSV_FIELD_NAMES = (
    "Reference Position",
    "Type",
    "Length",
    "Reference",
    "Allele",
    "Linkage",
    "Count",
    "Coverage",
    "Frequency",
    "Forward/reverse balance",
    "Average quality",
    "Overlapping annotations",
    "Coding region change",
    "Amino acid change",
)


class _FileArtifacts:

    __slots__ = ["file_path", "workding_id", "exec_seq", "cwd", "res_dir"]

    def __init__(
        self,
        *file_path: str,
        working_id: str = "",
        exec_seq: str = "",
        cwd: str = None,
        res_dir: str = None
    ) -> None:
        self.file_path: tuple[str] = file_path
        self.workding_id: str = working_id
        self.exec_seq: str = exec_seq
        self.cwd: str = cwd
        if res_dir is None:
            if cwd is None:
                raise ValueError("cwd and res_dir cannot both be None")
            res_dir = cwd
        self.res_dir: str = res_dir if res_dir is not None else cwd

    def __str__(self) -> str:
        return f"file_path: {self.file_path}, cwd: {self.cwd}"

    def __repr__(self) -> str:
        return str(self)

    def __len__(self) -> int:
        return len(self.file_path)

    def path_from_cwd(self) -> tuple:
        res = tuple(fp if os.path.isabs(fp)
                    else os.path.relpath(fp, self.cwd)
                    for fp in self.file_path)
        return res

    def exist(self) -> bool:
        if len(self.file_path) == 0:
            return False
        return all(os.path.exists(fp) for fp in self.file_path)

    def file_prefix(self) -> str:
        return os.path.join(self.cwd, self.workding_id)

    def relative_path(self, target_path: str) -> str:
        return os.path.relpath(target_path, self.cwd)

    def coupled_files(self, *exts: str, exec_name: str) -> "_FileArtifacts":
        middle = self.exec_seq
        middle += f"_{exec_name.replace(' ', '_')}" if exec_name else ""
        return _FileArtifacts(
            *tuple(f"{self.file_prefix()}{middle}{ext}" for ext in exts),
            working_id=self.workding_id,
            exec_seq=middle,
            cwd=self.cwd,
            res_dir=self.res_dir
        )

    def create_cwd(self) -> None:
        if os.path.exists(self.cwd):
            print(f"Delete existing {self.cwd}")
            shutil.rmtree(self.cwd)
        os.makedirs(self.cwd)

    def result_exists(self) -> bool:
        return os.path.exists(self.result_file())

    def result_file(self) -> None:
        return os.path.join(self.res_dir, f"{self.workding_id}.csv")

    def _write_result(self, workflow_res: tuple) -> None:
        with open(self.result_file(), "w", newline="") as f:
            csv_writer = csv.DictWriter(f, fieldnames=CSV_FIELD_NAMES)
            csv_writer.writeheader()
            csv_writer.writerows(workflow_res)


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
        if not self.input_files.exist():
            with open(self._log_file(), "a") as f:
                f.write(f"Input files for {str(self)} not found:\n")
                f.write(f"{', '.join(self.input_files.file_path)}")
            return _FileArtifacts(
                cwd=self.input_files.cwd,
                res_dir=self.input_files.res_dir
            )
        print("\n", str(self))
        try:
            with subprocess.Popen(
                self.cmd,
                text=True,
                cwd=self.input_files.cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            ) as p:
                self.stdout, self.stderr = p.communicate()
        except subprocess.SubprocessError as e:
            with open(self._log_file(), "a") as f:
                f.write(f"{str(self)} failed: {e}\n")
                f.write(f"Working directory: {self.input_files.cwd}\n")
                f.write(f"Executed command: {str(self)}\n")
        if p.returncode:
            self.append_log()
        else:
            self._post_execution()
        return self.output_files

    def append_log(self) -> None:
        with open(self._log_file(), "a") as f:
            f.write(f"{str(self)}\n")
            f.write(f"{self.exec_name} stdout:\n{self.stdout}\n")
            f.write(f"{self.exec_name} stderr:\n{self.stderr}\n")

    @abstractmethod
    def _post_execution(self) -> None:
        return NotImplemented

    # @abstractmethod
    # def _output_files(self) -> _FileArtifacts:
    #     return NotImplemented

    def _log_file(self) -> str:
        res = f"{self.input_files.file_prefix()}_log.txt"
        return res
