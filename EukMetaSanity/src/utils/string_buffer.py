import ctypes


class StringBuffer:
    def __init__(self, output_path=None, buf_size=4e4, initial=""):
        assert buf_size > 0
        assert isinstance(initial, str)
        assert output_path is None or isinstance(output_path, str)
        # Initialize buffer and position in buffer
        self._buf_size = int(float(buf_size)) - 1
        self._data = ctypes.create_string_buffer(b"\0" * (self._buf_size + 1))
        self._pos = 0
        # Open output buffer
        self._output = None
        if output_path is not None:
            self._output = open(output_path, "w")
        # Add initial value, if present
        self._add_to_buffer(initial)

    def _reinitialize_buffer(self):
        self._pos = 0

    def _to_str(self) -> str:
        return "".join([chr(_v) for _v in self._data[:self._pos]])

    def _add_to_buffer(self, _value: str):
        _len_val = len(_value)
        assert _len_val < self._buf_size - 1
        if self._pos + len(_value) >= self._buf_size - 1:
            self._flush_buffer()
        for i in range(_len_val):
            self._data[self._pos + i] = ord(_value[i])
        self._pos += _len_val

    def _flush_buffer(self, reinit=True):
        if self._output is not None:
            self._output.write(self._to_str())
        if reinit:
            self._reinitialize_buffer()

    def __iadd__(self, other):
        self._add_to_buffer(other)
        return self

    def __repr__(self) -> str:
        return self._to_str()

    def add(self, value: str):
        self._add_to_buffer(value)

    def write(self, endline: str = None):
        self._flush_buffer()
        if endline is not None:
            self._output.write(endline)

    def __str__(self) -> str:
        return self.__repr__()

    def __len__(self) -> int:
        return self._pos
