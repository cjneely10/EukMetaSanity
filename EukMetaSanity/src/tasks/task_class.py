from abc import ABC, abstractmethod


class Task(ABC):
    def __init__(self, input_paths_dict, cfg, pm, record_id, db_name):
        # Store passed input flag:input_path dict
        self._input_paths = input_paths_dict
        # Instantiate output dict variable
        self.output_paths = {}
        # Store threads and workers
        self._threads_pw = int(cfg.config.get(db_name, {}).get("THREADS", 1)),
        self._workers = int(cfg.config.get(db_name, {}).get("WORKERS", 1)),
        # Store path manager
        self._pm = pm
        # Store config manager
        self._cfg = cfg
        # Add name of db
        pm.add_dir(db_name)
        # Store working directory
        self._wdir = pm.get_dir(record_id, [db_name])
        super().__init__()

    @property
    def cfg(self):
        return self._cfg

    @property
    def wdir(self):
        return self._wdir

    @property
    def threads(self):
        return self._threads_pw

    @property
    def workers(self):
        return self._workers

    @property
    def pm(self):
        return self._pm

    @abstractmethod
    def run(self):
        assert self.output_paths is not None

    @abstractmethod
    def results(self):
        if self.output_paths is None:
            self.run()
        return self.output_paths


if __name__ == "__main__":
    pass
