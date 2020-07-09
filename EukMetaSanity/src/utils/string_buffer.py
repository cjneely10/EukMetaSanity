import numpy as np


class Buffer:
    def __init__(self, output: str = "/dev/stdout", buf_size: int = int(1e6), initial: str = ""):
        # Initialize buffer and position in buffer
        self._initialize_buffer(buf_size)
        self._buf_size = buf_size
        # Open output buffer
        self._output = open(output, "w")
        # Add initial value, if present
        self._add_to_buffer(initial)

    def _initialize_buffer(self, buf_size: int):
        self._data = np.array([0] * buf_size, dtype="int8")
        self._pos = 0

    def _to_str(self) -> str:
        return "".join(chr(_char) for _char in self._data[:self._pos])

    @property
    def buffer(self) -> str:
        return self._to_str()

    def _add_to_buffer(self, _value: str):
        if self._pos + len(_value) >= self._buf_size:
            self._flush_buffer()

        for i, char in enumerate(_value):
            self._data[self._pos + i] = ord(char)
        self._pos += len(_value)

    def _flush_buffer(self):
        self._output.write(self._to_str())
        self._initialize_buffer(self._buf_size)

    def __iadd__(self, other):
        self._add_to_buffer(other)
        return self

    def __repr__(self):
        return self._to_str()
