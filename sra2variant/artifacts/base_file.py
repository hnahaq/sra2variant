import os
import csv
import shutil
import logging


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

    def create_cwd(self) -> bool:
        if os.path.exists(self.result_file()):
            logging.warning(
                f"{self.workding_id} result already existed in {self.result_file()}")
            return False
        if os.path.exists(self.cwd):
            logging.warning(f"{self.workding_id} delete existing {self.cwd}")
            shutil.rmtree(self.cwd)
        os.makedirs(self.cwd)
        return True

    def result_file(self) -> None:
        return os.path.join(self.res_dir, f"{self.workding_id}.csv")

    def log_file(self) -> str:
        return f"{self.file_prefix()}_log.txt"

    def _write_result(self, workflow_res: tuple) -> None:
        with open(self.result_file(), "w", newline="") as f:
            csv_writer = csv.DictWriter(f, fieldnames=CSV_FIELD_NAMES)
            csv_writer.writeheader()
            csv_writer.writerows(workflow_res)
